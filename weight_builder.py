from pathlib import Path
import argparse
import json
import calendar

import pandas as pd


def is_seoul_sig(sig):
    s = str(sig).strip()
    return s.startswith("11")


# Map Korean weekday to index, with Monday = 0, ..., Sunday = 6
YOIL_TO_INDEX = {
    "월": 0,  # Monday
    "화": 1,  # Tuesday
    "수": 2,  # Wednesday
    "목": 3,  # Thursday
    "금": 4,  # Friday
    "토": 5,  # Saturday
    "일": 6,  # Sunday
}


def build_hourly_weights_seoul(
    data_dir,
    pattern="data_202003_*h.csv",
):
    #OLD VERSION: returns weights_by_hour[0..23].
    data_path = Path(data_dir)
    weights_by_hour: dict[int, dict[int, dict[int, float]]] = {}

    for path in sorted(data_path.glob(pattern)):
        hour_str = path.stem.split("_")[-1].replace("h", "")
        hour = int(hour_str)

        df = pd.read_csv(path)

        mask_seoul = (
            df["출발 시군구 코드"].astype(str).str.startswith("11")
            & df["도착 시군구 코드"].astype(str).str.startswith("11")
        )
        df = df[mask_seoul].copy()

        if df.empty:
            weights_by_hour[hour] = {}
            continue

        move = pd.to_numeric(
            df["이동인구(합)"].astype(str).str.replace("*", "", regex=False),
            errors="coerce",
        ).fillna(0.0)

        df = df.assign(이동인구_num=move)

        # Aggregate per OD pair (summing over weekday, gender, age, type)
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
            # Average per day in month (31 days in March)
            hour_weights[sig_out][sig_in] = w / 31.0

        weights_by_hour[hour] = hour_weights

    return weights_by_hour


def build_hourly_weekday_weights_seoul(
    data_dir,
    pattern="data_202003_*h.csv",
):
    #NEW VERSION: returns weights_by_hour_of_week[0..167].
    
    data_path = Path(data_dir)

    # Pre-allocate all 168 hours as empty dicts so keys 0..167 always exist
    weights_by_hour_of_week: dict[int, dict[int, dict[int, float]]] = {
        h: {} for h in range(168)
    }

    year = None
    month = None
    weekday_counts = None  

    for path in sorted(data_path.glob(pattern)):
        # Extract hour from filename, e.g. data_202003_06h.csv -> 6
        hour_str = path.stem.split("_")[-1].replace("h", "")
        hour = int(hour_str)

        df = pd.read_csv(path)

        # Determine year/month & weekday_counts once
        if year is None or month is None:
            if "대상연월" not in df.columns:
                raise ValueError("Column '대상연월' not found in CSV.")
            ym = int(df["대상연월"].iloc[0])
            year = ym // 100
            month = ym % 100

            # Count how many times each weekday occurs in that month
            cal = calendar.monthcalendar(year, month)
            weekday_counts = {i: 0 for i in range(7)}
            for week in cal:
                for i, day in enumerate(week):  # i: 0=Mon .. 6=Sun
                    if day != 0:
                        weekday_counts[i] += 1

        # Keep only rows where both origin & destination are Seoul SIGs
        mask_seoul = (
            df["출발 시군구 코드"].astype(str).str.startswith("11")
            & df["도착 시군구 코드"].astype(str).str.startswith("11")
        )
        df = df[mask_seoul].copy()

        if df.empty:
            # Nothing to add for this hour; keys remain empty dicts
            continue

        move = pd.to_numeric(
            df["이동인구(합)"].astype(str).str.replace("*", "", regex=False),
            errors="coerce",
        ).fillna(0.0)

        df = df.assign(이동인구_num=move)

        grouped = (
            df.groupby(["요일", "출발 시군구 코드", "도착 시군구 코드"])["이동인구_num"]
              .sum()
        )

        for (yoil, sig_out, sig_in), w in grouped.items():

            weekday_index = YOIL_TO_INDEX[yoil]
            hour_of_week = weekday_index * 24 + hour

            sig_out = int(sig_out)
            sig_in = int(sig_in)
            w = float(w)

            # Average per occurrence of that weekday in the month
            count_days = weekday_counts[weekday_index]

            #per_day_weight = w / float(count_days)
            per_day_weight = w
            
            if sig_out not in weights_by_hour_of_week[hour_of_week]:
                weights_by_hour_of_week[hour_of_week][sig_out] = {}

            prev = weights_by_hour_of_week[hour_of_week][sig_out].get(sig_in, 0.0)
            weights_by_hour_of_week[hour_of_week][sig_out][sig_in] = prev + per_day_weight

    return weights_by_hour_of_week


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
        help=(
            "Optional path to save weights as JSON. "
            "Keys: hour_of_week (0..167) -> sig_out -> sig_in -> weight."
        ),
    )

    args = parser.parse_args()

    weights = build_hourly_weekday_weights_seoul(args.data_dir)

    print("Hours-of-week processed:", sorted(k for k, v in weights.items() if v))

    if args.output:
        json_ready = {
            str(h): {
                str(o): {str(d): w for d, w in dests.items()}
                for o, dests in weights[h].items()
            }
            for h in sorted(weights.keys())
        }
        with open(args.output, "w", encoding="utf-8") as f:
            json.dump(json_ready, f, ensure_ascii=False, indent=2)
        print(f"\nSaved weights to {args.output}")


if __name__ == "__main__":
    main()

