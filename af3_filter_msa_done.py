#!/usr/bin/env python3
"""Migrate completed AlphaFast markers up and filter inputs before a (re)run.

For each ``<INPUT_DIR>/<name>.json`` input:

1. Any prior marker at
   ``<OUTPUT_DIR>/msa_output/partition_*/<name>/<name>_data.json`` is moved up
   to ``<OUTPUT_DIR>/<name>_data.json`` (the now-empty per-name dir is removed).
2. Inputs whose marker exists at ``<OUTPUT_DIR>/<name>_data.json`` are treated
   as done and skipped. The remaining inputs are symlinked into
   ``<FILTERED_DIR>`` so AlphaFast can be pointed at that directory and process
   only the unfinished cases.

Exit codes:
    0  — inputs remain; filtered_dir has been populated with symlinks
    10 — all inputs already complete; filtered_dir left empty
    other — error (bad args, input dir missing, etc.)
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys


def iter_input_names(input_dir):
    with os.scandir(input_dir) as it:
        for entry in it:
            name = entry.name
            if name.endswith(".json") and entry.is_file(follow_symlinks=True):
                yield name[:-5]


def reset_symlink_dir(path):
    """Empty a directory of symlinks (or create it). Refuses to recurse."""
    if os.path.isdir(path):
        with os.scandir(path) as it:
            for entry in it:
                if entry.is_symlink() or entry.is_file(follow_symlinks=False):
                    os.unlink(entry.path)
                else:
                    raise RuntimeError(
                        f"refusing to clean {path}: contains non-symlink entry "
                        f"{entry.name}"
                    )
    else:
        os.makedirs(path, exist_ok=True)


def migrate_markers(output_dir):
    """Move msa_output/partition_*/<name>/<name>_data.json -> <output_dir>/<name>_data.json."""
    part_root = os.path.join(output_dir, "msa_output")
    if not os.path.isdir(part_root):
        return 0
    moved = 0
    with os.scandir(part_root) as parts:
        for part in parts:
            if not (part.name.startswith("partition_") and part.is_dir(follow_symlinks=False)):
                continue
            with os.scandir(part.path) as it:
                for entry in it:
                    if not entry.is_dir(follow_symlinks=False):
                        continue
                    src = os.path.join(entry.path, f"{entry.name}_data.json")
                    if not os.path.isfile(src):
                        continue
                    dst = os.path.join(output_dir, f"{entry.name}_data.json")
                    shutil.move(src, dst)
                    moved += 1
                    try:
                        os.rmdir(entry.path)
                    except OSError:
                        pass
    return moved


def collect_done_names(output_dir):
    """Names whose marker exists at <output_dir>/<name>_data.json."""
    done = set()
    if not os.path.isdir(output_dir):
        return done
    suffix = "_data.json"
    with os.scandir(output_dir) as it:
        for entry in it:
            if entry.name.endswith(suffix) and entry.is_file(follow_symlinks=False):
                done.add(entry.name[: -len(suffix)])
    return done


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("input_dir")
    ap.add_argument("output_dir")
    ap.add_argument("filtered_dir")
    args = ap.parse_args()

    input_dir = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output_dir)
    filtered_dir = os.path.abspath(args.filtered_dir)

    if not os.path.isdir(input_dir):
        print(f"[filter] input dir not found: {input_dir}", file=sys.stderr)
        return 2

    names = list(iter_input_names(input_dir))
    total = len(names)
    if total == 0:
        print(f"[filter] no *.json inputs in {input_dir}", file=sys.stderr)
        return 2

    moved = migrate_markers(output_dir)
    if moved:
        print(f"[filter] migrated {moved} marker(s) up to {output_dir}")

    done_names = collect_done_names(output_dir)
    todo = [n for n in names if n not in done_names]
    done = total - len(todo)
    print(f"[filter] total={total} done={done} todo={len(todo)}")

    reset_symlink_dir(filtered_dir)

    if not todo:
        return 10

    for name in todo:
        os.symlink(
            os.path.join(input_dir, f"{name}.json"),
            os.path.join(filtered_dir, f"{name}.json"),
        )

    print(f"[filter] symlinked {len(todo)} unfinished input(s) -> {filtered_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
