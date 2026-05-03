#!/usr/bin/env python3
"""
gem_hycal_matching.py — Python counterpart of analysis/scripts/gem_hycal_matching.C

Same pipeline (HyCal reco → GEM reco → straight-line target-vertex matching),
same trigger filter (trigger_bits == 0x100), same multi-file discovery
(glob / directory / single).  Difference vs. the ROOT script:

  * No ROOT.  Output is a flat per-match TSV / CSV table — one row per
    (event, HyCal cluster, GEM detector) tuple, with `det_id` distinguishing
    which GEM (0..3) won the match.

  * Best-match rule (HyCal cluster as baseline):
      For each (HC cluster, GEM detector) pair, keep at most ONE row —
      the GEM hit with the smallest 2D residual that's still inside the
      `--match-nsigma · σ_total` window.  A given GEM hit can win against
      multiple HC clusters (no GEM-side exclusivity).

      Same rule as the patched C++ script — both produce the same matches
      bit-for-bit when run on the same EVIO with the same parameters.

Matching geometry (lab frame, target at origin, beam along +z):

    σ_hc_face = sqrt((A/sqrt(E_GeV))^2 + (B/E_GeV)^2 + C^2)  [mm at HyCal face]
    σ_hc@gem  = σ_hc_face · (z_gem / z_hc)
    σ_gem     = gem_pos_res[det_id] mm                       (per detector)
    σ_total   = sqrt(σ_hc@gem² + σ_gem²)
    cut       = nsigma · σ_total

    (A, B, C) and gem_pos_res come from reconstruction_config.json:matching.

Output columns (one row per matched pair):

  event_num, trigger_bits,
  hc_idx, hc_x, hc_y, hc_z, hc_energy, hc_center, hc_nblocks, hc_sigma,
  det_id,
  gem_x, gem_y, gem_z,                 # lab/target-centered (mm)
  gem_x_local, gem_y_local,            # detector-frame (mm)
  gem_x_charge, gem_y_charge,          # X/Y cluster total ADC
  gem_x_peak,   gem_y_peak,            # X/Y cluster max-strip ADC
  gem_x_max_tb, gem_y_max_tb,          # time sample of max-ADC strip (int)
  gem_x_size,   gem_y_size,            # X/Y cluster strip count
  proj_x, proj_y, residual, sigma_total

Coordinates labelled "lab" are target-centered (mm); hc_z includes
shower-depth.  The "_local" coords are the GEM detector frame (no
rotation/translation), useful for per-detector hit maps.  Convert
gem_*_max_tb to ns by multiplying by the cluster config's ts_period
(default 25 ns).

Usage
-----
  # full run (glob — warns about any missing split):
  python analysis/pyscripts/gem_hycal_matching.py \\
      /data/stage6/prad_023867/prad_023867.evio.* match_023867.tsv

  # single split (debugging):
  python analysis/pyscripts/gem_hycal_matching.py \\
      /data/stage6/prad_023867/prad_023867.evio.00000 match_023867_seg0.tsv

  # CSV output, tighter cut, cap at 50k events:
  python analysis/pyscripts/gem_hycal_matching.py input.evio.* out.csv \\
      --csv --match-nsigma 2.0 --max-events 50000
"""

from __future__ import annotations

import argparse
import math
import sys
import time

# _common imports prad2py and prints a friendly error if it's missing —
# import it first so we don't surface the bare ImportError instead.
import _common as C
from prad2py import dec, det  # noqa: E402  (after _common, intentionally)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    C.add_common_args(ap)
    ap.add_argument("--match-nsigma", type=float, default=3.0,
                    help="Matching window in σ_total (default 3.0).")
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
    print(f"[setup] Match cut  : {args.match_nsigma:.2f} · sigma_total",
          flush=True)

    (pr_A, pr_B, pr_C), gem_pos_res, _ = C.load_matching_config()
    print(f"[setup] HC sigma(E)= sqrt(({pr_A:.3f}/sqrt(E_GeV))^2"
          f"+({pr_B:.3f}/E_GeV)^2+{pr_C:.3f}^2) mm", flush=True)
    print(f"[setup] GEM sigma  : {gem_pos_res} mm", flush=True)

    cols = [
        "event_num", "trigger_bits",
        "hc_idx", "hc_x", "hc_y", "hc_z", "hc_energy",
        "hc_center", "hc_nblocks", "hc_sigma",
        "det_id",
        "gem_x", "gem_y", "gem_z",
        "gem_x_local", "gem_y_local",
        "gem_x_charge", "gem_y_charge",
        "gem_x_peak",   "gem_y_peak",
        "gem_x_max_tb", "gem_y_max_tb",
        "gem_x_size",   "gem_y_size",
        "proj_x", "proj_y", "residual", "sigma_total",
    ]
    fh, write_row = C.open_table_writer(args.out_path, args.csv)
    if not args.no_header:
        write_row(cols)

    ch = dec.EvChannel()
    ch.set_config(p.cfg)

    t0 = time.monotonic()
    n_read = n_phys = n_kept = n_match = 0
    n_files_open = 0
    total_hc = 0
    total_gem = 0
    gem_per_det = [0, 0, 0, 0]

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

                    # Trigger filter — gate against n_phys (raw count) so
                    # max_events cuts behave identically to the C++ script.
                    if fadc_evt.info.trigger_bits != C.PHYSICS_TRIGGER_BITS:
                        if args.max_events > 0 and n_phys >= args.max_events:
                            done = True; break
                        continue
                    n_kept += 1

                    event_num    = int(fadc_evt.info.event_number)
                    trigger_bits = int(fadc_evt.info.trigger_bits)

                    # ---- HyCal: waveform → energy → cluster --------------
                    hc_raw = C.reconstruct_hycal(p, fadc_evt)

                    # Lab-frame HyCal hits with shower-depth applied to z
                    # (det.shower_depth is the prad2det helper bound from
                    # fdec::shower_depth — same calc as Replay.cpp).
                    hc_lab: list[tuple[float, float, float, float, int, int]] = []
                    for h in hc_raw:
                        z_local = det.shower_depth(h.center_id, h.energy)
                        x, y, z = p.hycal_xform.to_lab(h.x, h.y, z_local)
                        hc_lab.append(
                            (x, y, z, float(h.energy),
                             int(h.center_id), int(h.nblocks))
                        )
                    total_hc += len(hc_lab)

                    # ---- GEM: pedestal → CM → ZS → 1D + 2D --------------
                    C.reconstruct_gem(p, ssp_evt)

                    # Per-detector lab-frame hit lists, plus the raw GEMHit
                    # for charge / size / peak / timing lookup at write time.
                    # Tuple layout (positional, frozen):
                    #   0: x_lab     1: y_lab     2: z_lab
                    #   3: x_local   4: y_local
                    #   5: x_charge  6: y_charge
                    #   7: x_peak    8: y_peak
                    #   9: x_max_tb 10: y_max_tb
                    #  11: x_size   12: y_size
                    gem_lab: list[list[tuple]] = [[], [], [], []]
                    n_dets = min(p.gem_sys.get_n_detectors(), 4)
                    for d in range(n_dets):
                        raw = p.gem_sys.get_hits(d)
                        gem_per_det[d] += len(raw)
                        total_gem      += len(raw)
                        xform = p.gem_xforms[d]
                        for g in raw:
                            x, y, z = xform.to_lab(g.x, g.y)
                            gem_lab[d].append((
                                x, y, z,
                                float(g.x), float(g.y),
                                float(g.x_charge), float(g.y_charge),
                                float(g.x_peak),   float(g.y_peak),
                                int(g.x_max_timebin), int(g.y_max_timebin),
                                int(g.x_size),     int(g.y_size),
                            ))

                    # ---- best match per HC × GEM detector ---------------
                    for k, (hx, hy, hz, he, hc_center, hc_nblocks) in enumerate(hc_lab):
                        if hz <= 0.0:
                            continue
                        # σ_HC at HyCal face — see _common.hycal_pos_resolution
                        sig_face = C.hycal_pos_resolution(pr_A, pr_B, pr_C, he)

                        for d in range(n_dets):
                            gl = gem_lab[d]
                            if not gl:
                                continue
                            z_gem = gl[0][2]
                            if z_gem <= 0.0:
                                continue
                            scale         = z_gem / hz
                            proj_x        = hx * scale
                            proj_y        = hy * scale
                            sig_hc_at_gem = sig_face * scale
                            sig_gem       = (gem_pos_res[d] if d < len(gem_pos_res)
                                             else 0.1)
                            sig_total     = math.sqrt(
                                sig_hc_at_gem * sig_hc_at_gem + sig_gem * sig_gem)
                            cut           = args.match_nsigma * sig_total

                            best_gi = -1
                            best_dr = cut
                            for gi, g in enumerate(gl):
                                dx = g[0] - proj_x
                                dy = g[1] - proj_y
                                dr = math.sqrt(dx * dx + dy * dy)
                                if dr <= best_dr:
                                    best_dr = dr
                                    best_gi = gi
                            if best_gi < 0:
                                continue

                            g = gl[best_gi]
                            write_row([
                                event_num, trigger_bits,
                                k,
                                f"{hx:.4f}", f"{hy:.4f}", f"{hz:.4f}",
                                f"{he:.4f}",
                                hc_center, hc_nblocks,
                                f"{sig_face:.4f}",
                                d,
                                f"{g[0]:.4f}", f"{g[1]:.4f}", f"{g[2]:.4f}",
                                f"{g[3]:.4f}", f"{g[4]:.4f}",
                                f"{g[5]:.4f}", f"{g[6]:.4f}",
                                f"{g[7]:.4f}", f"{g[8]:.4f}",
                                g[9], g[10],
                                g[11], g[12],
                                f"{proj_x:.4f}", f"{proj_y:.4f}",
                                f"{best_dr:.4f}", f"{sig_total:.4f}",
                            ])
                            n_match += 1

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
    print("--- summary ---", flush=True)
    print(f"  EVIO files opened     : {n_files_open} / {len(p.evio_files)}")
    print(f"  EVIO records          : {n_read}")
    print(f"  physics events        : {n_phys}")
    print(f"  passed trig cut 0x100 : {n_kept}")
    print(f"  total HyCal clusters  : {total_hc}")
    print(f"  total GEM 2D hits     : {total_gem}  "
          f"(det0={gem_per_det[0]} det1={gem_per_det[1]} "
          f"det2={gem_per_det[2]} det3={gem_per_det[3]})")
    print(f"  matched rows written  : {n_match}")
    print(f"  elapsed (s)           : {elapsed:.2f}")
    print(f"  wrote                 : {args.out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
