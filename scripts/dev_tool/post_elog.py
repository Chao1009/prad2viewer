#!/usr/bin/env python3
"""
post_elog.py — manually post a saved local auto-report to the JLab logbook.

Mirrors the prad2_server logic:

  * Dedup check: GET <elog_url>/api/elog/entries?book=<book>&title=...&
      field=lognumber&field=title&limit=20      (server: handleElogCheck)
    Substring server-side, EXACT title match client-side.

  * Post:        PUT <elog_url>/incoming/prad2_report.xml with the XML
                 body, cert + key auth                 (server: handleElogPost)

Stdlib only — no `requests`, no `pip install`.

Usage:
  ./post_elog.py 24364                       # latest XML in run_024364/
  ./post_elog.py 024364                      # leading zeros OK
  ./post_elog.py path/to/report.xml          # explicit path
  ./post_elog.py 24364 --check-only          # only print dedup result
  ./post_elog.py 24364 --force               # skip dedup check
"""

import argparse
import json
import re
import ssl
import sys
import urllib.parse
import urllib.request
from pathlib import Path

# Defaults match clonfarm11 / clonpc19 setup. Override via flags or by
# pointing --config at a monitor_config.json.
DEFAULT_REPORTS_DIR = "/home/clasrun/prad2_daq/monitor/reports"
DEFAULT_CERT        = "/home/clasrun/prad2_daq/monitor/keys/elog-cert.pem"
DEFAULT_KEY         = "/home/clasrun/prad2_daq/monitor/keys/elog-key.pem"
DEFAULT_URL         = "https://logbooks.jlab.org"
DEFAULT_BOOK        = "PRADLOG"
TITLE_BASE          = "PRad2 Event Monitor Auto Report"
TIMEOUT_S           = 30


def auto_title(run: int) -> str:
    return f"Run #{run}: {TITLE_BASE}"


def find_xml_for_run(reports_dir: str, run: int):
    rundir = Path(reports_dir) / f"run_{run:06d}"
    if not rundir.is_dir():
        return None
    xmls = sorted(rundir.glob("report_*.xml"),
                  key=lambda p: p.stat().st_mtime, reverse=True)
    return xmls[0] if xmls else None


def run_from_path(path) -> int | None:
    m = re.search(r"run_0*(\d+)", str(path))
    return int(m.group(1)) if m else None


def make_ssl_ctx(cert: str, key: str) -> ssl.SSLContext:
    ctx = ssl.create_default_context()
    ctx.load_cert_chain(certfile=cert, keyfile=key)
    return ctx


def check_elog(url: str, book: str, cert: str, key: str, title: str):
    """Mirror handleElogCheck — exact-title match against the live logbook."""
    qs = urllib.parse.urlencode([
        ("book",  book),
        ("title", title),
        ("field", "lognumber"),
        ("field", "title"),
        ("limit", "20"),
    ])
    full = f"{url}/api/elog/entries?{qs}"
    ctx = make_ssl_ctx(cert, key)
    req = urllib.request.Request(full, method="GET")
    try:
        with urllib.request.urlopen(req, context=ctx, timeout=TIMEOUT_S) as r:
            body = r.read().decode("utf-8", errors="replace")
            code = r.status
    except urllib.error.HTTPError as e:
        return {"checked": False, "detail": f"HTTP {e.code}"}
    except Exception as e:
        return {"checked": False, "detail": f"network: {e}"}
    if code != 200:
        return {"checked": False, "detail": f"HTTP {code}"}
    try:
        j = json.loads(body)
    except json.JSONDecodeError as e:
        return {"checked": False, "detail": f"json: {e}"}
    if j.get("stat") != "ok":
        return {"checked": False,
                "detail": j.get("message", "non-ok stat")}
    entries = j.get("data", {}).get("entries", [])
    for e in entries:
        if e.get("title") == title:
            return {"checked": True, "exists": True,
                    "lognumber": e.get("lognumber"),
                    "matched_count": len(entries)}
    return {"checked": True, "exists": False,
            "matched_count": len(entries)}


def post_elog(url: str, cert: str, key: str, xml_path: Path):
    """Mirror handleElogPost — HTTP PUT to /incoming/prad2_report.xml."""
    body = xml_path.read_bytes()
    ctx  = make_ssl_ctx(cert, key)
    full = f"{url}/incoming/prad2_report.xml"
    req  = urllib.request.Request(full, data=body, method="PUT",
                                  headers={"Content-Type": "application/xml"})
    try:
        with urllib.request.urlopen(req, context=ctx, timeout=TIMEOUT_S) as r:
            return r.status, r.read().decode("utf-8", errors="replace")
    except urllib.error.HTTPError as e:
        return e.code, (e.read().decode("utf-8", errors="replace") if e.fp else "")
    except Exception as e:
        return 0, f"network: {e}"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("target",
        help="run number (e.g. 24364 or 024364) or path to a report.xml")
    ap.add_argument("--reports-dir", default=DEFAULT_REPORTS_DIR)
    ap.add_argument("--cert", default=DEFAULT_CERT)
    ap.add_argument("--key",  default=DEFAULT_KEY)
    ap.add_argument("--url",  default=DEFAULT_URL)
    ap.add_argument("--book", default=DEFAULT_BOOK)
    ap.add_argument("--run",  type=int,
        help="override run number when target is an XML path")
    ap.add_argument("--force", action="store_true",
        help="skip the dedup check (post even if already in elog)")
    ap.add_argument("--check-only", action="store_true",
        help="print dedup result and exit; do not post")
    args = ap.parse_args(argv)

    # Resolve XML path + run number
    if Path(args.target).is_file():
        xml_path = Path(args.target)
        run = args.run if args.run is not None else run_from_path(xml_path)
    else:
        try:
            run = args.run if args.run is not None else int(args.target.lstrip("0") or "0")
        except ValueError:
            sys.exit(f"could not parse run number from '{args.target}'")
        xml_path = find_xml_for_run(args.reports_dir, run)
        if xml_path is None:
            sys.exit(f"no XML found in {args.reports_dir}/run_{run:06d}")
    if not run:
        sys.exit("could not determine run number")

    title = auto_title(run)
    print(f"title : {title}")
    print(f"xml   : {xml_path}")

    if not args.force:
        print("checking elog...")
        c = check_elog(args.url, args.book, args.cert, args.key, title)
        if not c.get("checked"):
            print(f"  dedup unavailable: {c.get('detail')}")
            if not args.check_only:
                # Match the server's fail-closed policy on uncertain dedup.
                sys.exit("dedup couldn't run — refusing to post (re-run with --force to override)")
        elif c.get("exists"):
            print(f"  ALREADY IN ELOG: lognumber={c.get('lognumber')}  (skipping)")
            sys.exit(0)
        else:
            print(f"  not in elog ({c.get('matched_count')} fuzzy hits)")
        if args.check_only:
            return

    print("posting...")
    code, body = post_elog(args.url, args.cert, args.key, xml_path)
    print(f"HTTP {code}")
    if body:
        print(body[:1000])
    sys.exit(0 if code in (200, 201) else 1)


if __name__ == "__main__":
    main()
