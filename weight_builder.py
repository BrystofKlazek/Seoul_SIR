#!/usr/bin/env python3
#!/usr/bin/env python3
"""
Build hourly OD weights BETWEEN SEOUL DISTRICTS ONLY.

Input files are expected to look like:
    data_202003_00h.csv
    data_202003_01h.csv
    ...
    data_202003_23h.csv

Each file must have at least these columns:
    - 출발 시군구 코드
    - 도착 시군구 코드
    - 이동인구(합)

Output (Python object):
    weights[hour][sig_out][sig_in] = total_flow

Where:
    - hour is an int 0..23
    - sig_out, sig_in are ints (시군구 코드)
    - total_flow is float (sum of 이동인구(합) over all rows)

Only pairs where BOTH sig_out and sig_in are Seoul codes (11xxx) are kept.
"""

from pathlib import Path
import argparse
import json

import pandas as pd


def is_seoul_sig(sig) -> bool:
    """
    Return True if this 시군구 코드 is a Seoul code.

    Seoul districts have codes starting with '11', e.g. 11010, 11110, ...
    """
    s = str(sig).strip()
    return s.startswith("11")


def build_hourly_weights_seoul(data_dir: str,
                               pattern: str = "data_202003_*h.csv") -> dict:
    """
    Scan data_dir for hourly CSVs and return a nested dict:

        weights[hour][sig_out][sig_in] = weight

    including ONLY trips where both origin and destination
    are Seoul SIGs (11xxx codes).
    """
    data_path = Path(data_dir)
    weights_by_hour: dict[int, dict[int, dict[int, float]]] = {}

    for path in sorted(data_path.glob(pattern)):
        # Extract hour from filename, e.g. data_202003_06h.csv -> 6
        hour_str = path.stem.split("_")[-1].replace("h", "")
        hour = int(hour_str)

        df = pd.read_csv(path)

        # Keep only rows where both origin & destination are Seoul SIGs
        mask_seoul = (
            df["출발 시군구 코드"].astype(str).str.startswith("11")
            & df["도착 시군구 코드"].astype(str).str.startswith("11")
        )
        df = df[mask_seoul].copy()

        if df.empty:
            # No Seoul↔Seoul trips in this hour
            weights_by_hour[hour] = {}
            continue

        # Clean 이동인구(합):
        #   - remove '*' etc.
        #   - convert to numeric (errors -> NaN -> 0.0)
        move = pd.to_numeric(
            df["이동인구(합)"].astype(str).str.replace("*", "", regex=False),
            errors="coerce",
        ).fillna(0.0)

        df = df.assign(이동인구_num=move)

        # Aggregate per OD pair
        grouped = (
            df.groupby(["출발 시군구 코드", "도착 시군구 코드"])["이동인구_num"]
              .sum()
        )

        hour_weights: dict[int, dict[int, float]] = {}
        for (sig_out, sig_in), w in grouped.items():
            sig_out = int(sig_out)
            sig_in = int(sig_in)
            w = float(w)

            if sig_out not in hour_weights:
                hour_weights[sig_out] = {}
            hour_weights[sig_out][sig_in] = w/31

        weights_by_hour[hour] = hour_weights

    return weights_by_hour


def main():
    parser = argparse.ArgumentParser(
        description="Build hourly Seoul-only OD weight dicts from CSV files."
    )
    parser.add_argument(
        "data_dir",
        help="Directory containing data_202003_??h.csv files",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Optional path to save weights as JSON "
             "(keys: hour -> sig_out -> sig_in -> weight).",
    )

    args = parser.parse_args()

    weights = build_hourly_weights_seoul(args.data_dir)

    # Tiny summary
    print("Hours processed:", sorted(weights.keys()))
    for h in sorted(weights.keys()):
        num_edges = sum(len(dests) for dests in weights[h].values())
        print(f"  hour {h:02d}: {num_edges} Seoul↔Seoul OD pairs")

    if args.output:
        # Convert hour keys to strings for JSON
        json_ready = {
            str(h): {
                str(o): {str(d): w for d, w in dests.items()}
                for o, dests in weights[h].items()
            }
            for h in weights
        }
        with open(args.output, "w", encoding="utf-8") as f:
            json.dump(json_ready, f, ensure_ascii=False, indent=2)
        print(f"\nSaved weights to {args.output}")


if __name__ == "__main__":
    main()

