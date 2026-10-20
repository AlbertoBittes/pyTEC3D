"""Pytest configuration: run the tests without writing bytecode into the source tree."""

import sys

sys.dont_write_bytecode = True


def pytest_terminal_summary(terminalreporter, exitstatus, config):
    """Table of the APDL benchmarks (tests/test_apdl_benchmarks.py): worst relative L2 error per case and field.

    ``err`` is the worst error over the asserted sets (``set``: where it occurs); ``info`` the worst error over the sets
    where the field is only reported (e.g. Jt on steps with theta < 1); a field never asserted has no ``err`` and the
    status ``info``. With ``PYTEC3D_APDL_REPORT=<file.json>`` the
    table is also written as JSON (with the reference set and the inputs) and as Markdown (``<file>.md``).
    """
    import os

    module = sys.modules.get("test_apdl_benchmarks")
    rows = getattr(module, "SUMMARY", None) if module else None
    if module is None:
        return
    info = getattr(module, "RUN_INFO", {})
    tr = terminalreporter
    if not rows:
        if info.get("skipped"):
            tr.write_sep("=", "APDL benchmarks skipped")
            tr.write_line(info["skipped"])
        return
    tr.write_sep("=", "APDL benchmark summary (relative L2 error pyTEC3D vs MAPDL)")
    tr.write_line(f"references {info.get('references')}, inputs {info.get('inputs_fingerprint')}")
    tr.write_line(f"{'case':30s} {'field':14s} {'err':>10s} {'rtol':>8s} {'set':>4s} {'info':>10s}  status")
    last = None
    for r in rows:
        case = r["case"] if r["case"] != last else ""
        last = r["case"]
        info_err = f"{r['info_err']:10.2e}" if r["info_err"] is not None else f"{'':>10s}"
        err = f"{r['err']:10.2e}" if r["err"] is not None else f"{'-':>10s}"
        tr.write_line(f"{case:30s} {r['field']:14s} {err} {r['rtol']:8.0e} {r['set'] or '':>4} "
                      f"{info_err}  {_status(r)}", red=not r["passed"], green=r["passed"])
    cases = sorted({r["case"] for r in rows})
    failed = sorted({r["case"] for r in rows if not r["passed"]})
    tr.write_line(f"{len(cases) - len(failed)} / {len(cases)} cases within tolerance"
                  + (f"; failing: {', '.join(failed)}" if failed else ""))
    path = os.environ.get("PYTEC3D_APDL_REPORT")
    if path:
        md = write_apdl_report(path, rows, info)
        tr.write_line(f"report written to {path} and {md}")


def _status(r):
    """ok / FAIL (with the error of a failed run) / info (field never asserted: reported only)."""
    if not r["passed"]:
        return "FAIL" + (f": {r['error']}" if r.get("error") else "")
    return "info" if r["err"] is None else "ok"


def _pytec3d_revision():
    import os
    import subprocess

    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    try:
        rev = subprocess.run(["git", "-C", root, "describe", "--always", "--dirty"], capture_output=True, text=True,
                             check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        rev = "unknown"
    return rev


def write_apdl_report(path, rows, info):
    """APDL benchmark report: ``path`` (JSON) and the same path with extension ``.md`` (Markdown). Returns the latter."""
    import datetime
    import json
    import os

    cases = sorted({r["case"] for r in rows})
    failed = sorted({r["case"] for r in rows if not r["passed"]})
    report = {"date": datetime.datetime.now().isoformat(timespec="seconds"), "pytec3d": _pytec3d_revision(),
              **{k: v for k, v in info.items() if k != "skipped"},
              "cases": len(cases), "failed": failed, "rows": rows}
    with open(path, "w") as f:
        json.dump(report, f, indent=1)

    lines = [
        "# APDL benchmark assessment", "",
        f"**{len(cases) - len(failed)} / {len(cases)} cases within tolerance**"
        + (f" (failing: {', '.join(failed)})" if failed else ""), "",
        "| | |", "|---|---|",
        f"| date | {report['date']} |",
        f"| pyTEC3D | `{report['pytec3d']}` |",
        f"| MAPDL references | `{info.get('references')}` ({info.get('solver')}) |",
        f"| inputs fingerprint | `{info.get('inputs_fingerprint')}` |",
        f"| model data | `{info.get('model_data')}` |",
        f"| time history | `{info.get('time_history')}` |", "",
        "Relative L2 error ||pyTEC3D - MAPDL|| / ||MAPDL||: `err` is the worst over the asserted result sets (`set`),",
        "`info` the worst over the sets where the field is only reported.", "",
        "| case | field | err | rtol | set | info | status |", "|---|---|---:|---:|---:|---:|---|",
    ]
    last = None
    for r in rows:
        case = f"`{r['case']}`" if r["case"] != last else ""
        last = r["case"]
        info_err = f"{r['info_err']:.2e}" if r["info_err"] is not None else ""
        err = f"{r['err']:.2e}" if r["err"] is not None else "-"
        status = _status(r).replace("FAIL", "**FAIL**", 1)
        lines.append(f"| {case} | {r['field']} | {err} | {r['rtol']:.0e} | {r['set'] or ''} | {info_err} | {status} |")
    md = os.path.splitext(path)[0] + ".md"
    with open(md, "w") as f:
        f.write("\n".join(lines) + "\n")
    return md
