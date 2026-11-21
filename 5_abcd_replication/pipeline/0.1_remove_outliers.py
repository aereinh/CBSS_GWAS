import argparse
import numpy as np

def main():
    parser = argparse.ArgumentParser(
        description="Replace outliers in fibersClust txt with NA (after first two columns)"
    )
    parser.add_argument("--in-txt", required=True, help="Input fibersClust_*.txt")
    parser.add_argument("--out-txt", required=True, help="Output cleaned txt")
    parser.add_argument("--threshold", type=float, required=True,
                        help="Absolute value threshold; values > threshold set to NA")
    args = parser.parse_args()

    thresh = args.threshold

    with open(args.in_txt, "r") as fin, open(args.out_txt, "w") as fout:
        header = fin.readline()
        fout.write(header)

        for line in fin:
            parts = line.rstrip().split()
            if len(parts) <= 2:
                fout.write(line)
                continue

            # first two columns untouched
            out = parts[:2]

            # numeric part
            for v in parts[2:]:
                try:
                    x = float(v)
                    if (not np.isfinite(x)) or abs(x) > thresh:
                        out.append("NA")
                    else:
                        out.append(v)
                except ValueError:
                    out.append("NA")

            fout.write(" ".join(out) + "\n")

if __name__ == "__main__":
    main()