"""
Write the provenance table from the configuration actually used.

Protocol Section 12.0 requires that software version, database release and
every parameter be recorded. Transcribing them by hand after the run is how
protocols drift from what was executed, so the table is generated from the
config file the workflow read.
"""

import os
import platform
import subprocess
import sys
from datetime import datetime, timezone

import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import write_tsv

FIELDS = ["category", "item", "value", "source"]


def git_rev():
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], stderr=subprocess.DEVNULL,
            text=True).strip()
    except Exception:
        return "not_a_git_checkout"


def flatten(prefix, obj, out, source):
    if isinstance(obj, dict):
        for k, v in obj.items():
            flatten(f"{prefix}.{k}" if prefix else str(k), v, out, source)
    elif isinstance(obj, (list, tuple)):
        out.append((prefix, ",".join(str(x) for x in obj), source))
    else:
        out.append((prefix, "null" if obj is None else str(obj), source))


def main(cfg_path, tool, out_path):
    with open(cfg_path) as fh:
        cfg = yaml.safe_load(fh)

    rows = [
        {"category": "run", "item": "timestamp_utc",
         "value": datetime.now(timezone.utc).isoformat(timespec="seconds"),
         "source": "runtime"},
        {"category": "run", "item": "workflow_git_revision",
         "value": git_rev(), "source": "git"},
        {"category": "run", "item": "host", "value": platform.node(),
         "source": "runtime"},
        {"category": "run", "item": "platform", "value": platform.platform(),
         "source": "runtime"},
        {"category": "run", "item": "python", "value": platform.python_version(),
         "source": "runtime"},
        {"category": "run", "item": "independent_predictor_used",
         "value": tool, "source": "config"},
    ]

    for section, label in (("versions", "software"),
                           ("databases", "database"),
                           ("genomad", "parameter"),
                           ("checkv", "parameter"),
                           ("concordance", "parameter"),
                           ("quality", "parameter"),
                           ("votu", "parameter"),
                           ("poppunk", "parameter"),
                           ("threads", "resource"),
                           ("pide", "parameter"),
                           ("phispy", "parameter")):
        if section not in cfg:
            continue
        flat = []
        flatten(section, cfg[section], flat, "config")
        for item, value, src in flat:
            rows.append({"category": label, "item": item, "value": value,
                         "source": src})

    # Record the database directory modification time as a proxy for release
    # date where the release string was not filled in. A proxy labelled as a
    # proxy is more use than a blank.
    for name, path in (cfg.get("databases") or {}).items():
        if path and os.path.exists(path):
            mt = datetime.fromtimestamp(os.path.getmtime(path),
                                        timezone.utc).date().isoformat()
            rows.append({"category": "database",
                         "item": f"databases.{name}.mtime_proxy",
                         "value": mt, "source": "filesystem"})
        else:
            rows.append({"category": "database",
                         "item": f"databases.{name}.present",
                         "value": "MISSING", "source": "filesystem"})

    write_tsv(out_path, rows, FIELDS)
    missing = [r for r in rows if r["value"] == "MISSING"]
    if missing:
        sys.stderr.write(
            f"WARNING {len(missing)} configured database paths do not exist "
            f"on this host\n")
    sys.stderr.write(f"provenance: {len(rows)} entries\n")


if __name__ == "__main__":
    try:
        sm = snakemake  # noqa: F821
    except NameError:
        sm = None
    if sm is not None:
        main(sm.params.cfg, sm.params.tool, sm.output.tsv)
    else:
        main(sys.argv[1], sys.argv[2], sys.argv[3])
