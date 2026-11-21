#!/usr/bin/env python3
"""Rewrite the 'id' dataset of a HEIG HDF5 file (images.h5 or ldr_loco_preds.h5)
with full-length FID/IID strings.

HEIG truncates the long ABCD IDs (NDAR_INVxxxxxxxx) when it writes HDF5 files. This copies
the correct IDs from the first two columns of a whitespace-delimited text file with a
header (fibersClust_<IDX>.txt or <TAG>_ldr_top<N>.txt) whose rows are in the same order.
Safe to run more than once.
"""
import argparse
import os
import sys

import h5py
import numpy as np
import pandas as pd


def as_str(x):
    return x.decode("utf-8") if isinstance(x, (bytes, bytearray)) else str(x)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--h5", required=True)
    ap.add_argument("--ids-txt", required=True)
    args = ap.parse_args()

    for fp in (args.h5, args.ids_txt):
        if not os.path.exists(fp):
            sys.exit(f"ERROR: missing {fp}")

    df = pd.read_csv(args.ids_txt, sep=r"\s+", header=0, usecols=[0, 1], dtype=str)
    ids = df.fillna("").astype(str).to_numpy()
    if np.any(ids == ""):
        sys.exit(f"ERROR: blank FID/IID in {args.ids_txt}")

    with h5py.File(args.h5, "r+") as f:
        if "id" in f:
            old = f["id"][:]
            if old.shape[0] != ids.shape[0]:
                sys.exit(f"ERROR: row mismatch: {args.h5} has {old.shape[0]} ids, "
                         f"{args.ids_txt} has {ids.shape[0]}")
            old_iid = np.array([as_str(r[-1]) if np.ndim(r) else as_str(r) for r in old])
            # truncated IDs should be prefixes of the full IDs if the row order matches
            n_bad = int(sum(not new.startswith(o) for new, o in zip(ids[:, 1], old_iid)))
            if n_bad:
                print(f"WARNING: {n_bad}/{len(ids)} stored IDs are not a prefix of the new IDs "
                      f"(e.g. {old_iid[:3].tolist()} vs {ids[:3, 1].tolist()}); check row order.")
            del f["id"]
        else:
            print(f"No 'id' dataset in {args.h5}; creating it.")
        f.create_dataset("id", data=ids.astype(object), dtype=h5py.string_dtype(encoding="utf-8"))

    with h5py.File(args.h5, "r") as f:
        first = [as_str(x) for x in f["id"][0]]
    print(f"Rewrote 'id' in {args.h5} ({len(ids)} rows); first = {first}")


if __name__ == "__main__":
    main()
