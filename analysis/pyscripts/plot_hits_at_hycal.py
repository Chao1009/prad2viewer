#!/usr/bin/env python3
"""
plot_hits_at_hycal.py — Python counterpart of analysis/scripts/plot_hits_at_hycal.C

Same pipeline (HyCal reco → GEM reco → coord transform → straight-line
projection of GEM hits to z = hycal_z), same trigger filter
(trigger_bits == 0x100), same multi-file discovery (glob / directory /
single).  Difference vs. the ROOT script:

  * No ROOT.  Output is a flat per-hit TSV / CSV table — one row per
    hit (HyCal cluster centroid OR GEM hit projected to HyCal surface),
    with `kind` and `det_id` columns so you can split by source.

  * NOT matched-filtered.  This script dumps every hit the reconstruction
    produced, mirroring what the ROOT plot script feeds into its TH2Fs.
    For best-matched hits only, use gem_hycal_matching.py.

Output columns (one row per hit):

  event_num, trigger_bits, kind, det_id, x, y, z, energy

  kind   : "hycal" or "gem"
  det_id : 0..3 for GEM, -1 for HyCal
  x, y   : lab-frame mm at z = hycal_z
  z      : hycal_z (always) — kept for downstream uniformity
  energy : MeV for HyCal clusters, NaN (printed as empty in TSV) for GEM

The HyCal hits are placed on z = hycal_z exactly (no shower-depth
correction) so the (x, y) in this table matches the C++ plot script's
TH2F bin coordinates.

Usage
-----
  # full run:
  python analysis/pyscripts/plot_hits_at_hycal.py \\
      /data/stage6/prad_023867/prad_023867.evio.* hits_023867.tsv

  # single split:
  python analysis/pyscripts/plot_hits_at_hycal.py \\
      /data/stage6/prad_023867/prad_023867.evio.00000 hits_023867_seg0.tsv

  # CSV, cap events:
  python analysis/pyscripts/plot_hits_at_hycal.py input.evio.* out.csv \\
      --csv --max-events 50000

Plotting
--------
The output is plain rows; load with pandas / numpy and bin yourself, e.g.:

  import pandas as pd, matplotlib.pyplot as plt
  df = pd.read_csv("hits_023867.tsv", sep="\\t")
  for kind, sub in df.groupby("kind"):
      plt.hist2d(sub.x, sub.y, bins=260, range=[[-650, 650], [-650, 650]])
      plt.title(kind); plt.show()
"""

from __future__ import annotations

import argparse
import sys
import time

# _common imports prad2py and prints a friendly error if it's missing —
# import it first so we don't surface the bare ImportError instead.
import _common as C
from prad2py import dec  # noqa: E402  (after _common, intentionally)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    C.add_common_args(ap)
    args = ap.parse_args(argv)

    p = C.setup_pipeline(
        evio_path     = args.evio_path,
        max_events    = args.max_events,
        run_num       = args.run_num,
        gem_ped_file  = args.gem_ped_file,
        gem_cm_file   = args.gem_cm_file,
        hc_calib_file = args.hc_calib_file,
        daq_config    = args.daq_config,
        gem_map_file  = args.gem_map_file,
        hc_map_file   = args.hc_map_file,
    )

    cols = ["event_num", "trigger_bits", "kind", "det_id",
            "x", "y", "z", "energy"]
    fh, write_row = C.open_table_writer(args.out_path, args.csv)
    if not args.no_header:
        write_row(cols)

    ch = dec.EvChannel()
    ch.set_config(p.cfg)

    t0 = time.monotonic()
    n_read = n_phys = n_kept = 0
    n_files_open = 0
    n_hc_rows = 0
    n_gem_rows = 0

    try:
        for fpath in p.evio_files:
            if ch.open_auto(fpath) != dec.Status.success:
                print(f"[WARN] skip (cannot open): {fpath}", flush=True)
                continue
            n_files_open += 1
            print(f"[file {n_files_open}/{len(p.evio_files)}] {fpath}",
                  flush=True)

            done = False
            while ch.read() == dec.Status.success:
                n_read += 1
                if not ch.scan():
                    continue
                if ch.get_event_type() != dec.EventType.Physics:
                    continue

                for i in range(ch.get_n_events()):
                    decoded = ch.decode_event(i, with_ssp=True)
                    if not decoded["ok"]:
                        continue
                    n_phys += 1
                    fadc_evt = decoded["event"]
                    ssp_evt  = decoded["ssp"]

                    if fadc_evt.info.trigger_bits != C.PHYSICS_TRIGGER_BITS:
                        if args.max_events > 0 and n_phys >= args.max_events:
                            done = True; break
                        continue
                    n_kept += 1

                    event_num    = int(fadc_evt.info.event_number)
                    trigger_bits = int(fadc_evt.info.trigger_bits)
                    z_hycal      = p.geo.hycal_z

                    # ---- HyCal: build hits at z=0 (no shower depth) so
                    # the lab-frame transform lands them exactly on the
                    # HyCal face — same convention as plot_hits_at_hycal.C.
                    hc_raw = C.reconstruct_hycal(p, fadc_evt)
                    for h in hc_raw:
                        x, y, z = p.hycal_xform.to_lab(h.x, h.y)
                        write_row([event_num, trigger_bits, "hycal", -1,
                                   f"{x:.4f}", f"{y:.4f}", f"{z:.4f}",
                                   f"{float(h.energy):.4f}"])
                        n_hc_rows += 1

                    # ---- GEM: per-detector lab transform → project to
                    # z = hycal_z via straight-line target-vertex projection.
                    C.reconstruct_gem(p, ssp_evt)
                    n_dets = min(p.gem_sys.get_n_detectors(), 4)
                    for d in range(n_dets):
                        xform = p.gem_xforms[d]
                        for g in p.gem_sys.get_hits(d):
                            x, y, z = xform.to_lab(g.x, g.y)
                            x, y, z = C.project_to_z(x, y, z, z_hycal)
                            write_row([event_num, trigger_bits, "gem", d,
                                       f"{x:.4f}", f"{y:.4f}", f"{z:.4f}",
                                       ""])
                            n_gem_rows += 1

                    if args.max_events > 0 and n_phys >= args.max_events:
                        done = True; break

                if done:
                    break
                if n_phys > 0 and n_phys % 5000 == 0:
                    print(f"[progress] {n_phys} physics events", flush=True)

            ch.close()
            if done:
                break
    finally:
        fh.close()

    elapsed = time.monotonic() - t0
    avg_hc  = f" (avg {n_hc_rows  / n_kept:.2f} / kept event)" if n_kept else ""
    avg_gem = f" (avg {n_gem_rows / n_kept:.2f} / kept event)" if n_kept else ""
    print("--- summary ---", flush=True)
    print(f"  EVIO files opened     : {n_files_open} / {len(p.evio_files)}")
    print(f"  EVIO records          : {n_read}")
    print(f"  physics events        : {n_phys}")
    print(f"  passed trig cut 0x100 : {n_kept}")
    print(f"  HyCal rows written    : {n_hc_rows}{avg_hc}")
    print(f"  GEM rows written      : {n_gem_rows}{avg_gem}")
    print(f"  elapsed (s)           : {elapsed:.2f}")
    print(f"  wrote                 : {args.out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
