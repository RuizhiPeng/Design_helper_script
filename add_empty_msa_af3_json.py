#!/usr/bin/env python3
"""Transform AF3 input JSON files into the _data.json format, overwriting in place.

Usage: python transform_af3_json.py <input_folder>
"""
import argparse
import json
import sys
from pathlib import Path


def transform(data: dict) -> dict:
    out = {
        "dialect": data.get("dialect", "alphafold3"),
        "version": 4,
        "name": data["name"],
        "sequences": [],
        "modelSeeds": data.get("modelSeeds", []),
        "bondedAtomPairs": data.get("bondedAtomPairs", None),
        "userCCD": data.get("userCCD", None),
    }

    for seq in data.get("sequences", []):
        if "protein" in seq:
            p = seq["protein"]
            new_protein = {
                "id": p["id"],
                "sequence": p["sequence"],
                "modifications": p.get("modifications", []),
                "unpairedMsa": "",
                "pairedMsa": "",
                "templates": p.get("templates", []),
            }
            out["sequences"].append({"protein": new_protein})
        else:
            out["sequences"].append(seq)

    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_folder", help="Folder containing JSON files to transform in place.")
    args = parser.parse_args()

    folder = Path(args.input_folder)
    if not folder.is_dir():
        sys.exit(f"Error: {folder} is not a directory")

    files = sorted(folder.glob("*.json"))
    if not files:
        print(f"No .json files found in {folder}")
        return

    for path in files:
        with open(path) as f:
            data = json.load(f)
        transformed = transform(data)
        with open(path, "w") as f:
            json.dump(transformed, f, indent=2)
        print(f"Transformed: {path}")


if __name__ == "__main__":
    main()
