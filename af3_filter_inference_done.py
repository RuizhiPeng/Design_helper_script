#!/usr/bin/env python3
"""Stage pending AF3 inference inputs via symlinks.

Usage: af3_filter_inference_done.py <msa_dir> <out_dir> <fold_input_dir>

<msa_dir> is a flat directory of ``<sample>_data.json`` files produced by
the AF3 data pipeline. <out_dir> holds per-sample result folders; a sample
is considered done when ``<out_dir>/<sample>/<sample>_model.cif`` exists.
Pending samples are symlinked into <fold_input_dir>. Prints the number of
pending samples to stdout.
"""
import os
import sys
from concurrent.futures import ThreadPoolExecutor

SUFFIX = ".json"
DATA_SUFFIX = "_data"

def main():
    if len(sys.argv) != 4:
        print(f"usage: {sys.argv[0]} <msa_dir> <out_dir> <fold_input_dir>",
              file=sys.stderr)
        return 2

    msa_dir, out_dir, fold_input_dir = (os.path.abspath(p) for p in sys.argv[1:4])
    os.makedirs(fold_input_dir, exist_ok=True)

    done = set()
    if os.path.isdir(out_dir):
        with os.scandir(out_dir) as it:
            for e in it:
                if not e.is_dir(follow_symlinks=False):
                    continue
                if os.path.isfile(os.path.join(e.path, f"{e.name}_model.cif")):
                    done.add(e.name)

    pending = []
    with os.scandir(msa_dir) as it:
        for e in it:
            if not e.name.endswith(SUFFIX):
                continue
            if not e.is_file(follow_symlinks=True):
                continue
            sample = e.name[:-len(SUFFIX)]
            if sample.endswith(DATA_SUFFIX):
                sample = sample[:-len(DATA_SUFFIX)]
            if sample in done:
                continue
            pending.append((e.path, os.path.join(fold_input_dir, e.name)))

    def link(pair):
        src, dst = pair
        try:
            os.symlink(src, dst)
        except FileExistsError:
            os.unlink(dst)
            os.symlink(src, dst)

    if pending:
        with ThreadPoolExecutor(max_workers=16) as ex:
            for _ in ex.map(link, pending):
                pass

    print(len(pending))
    return 0


if __name__ == "__main__":
    sys.exit(main())
