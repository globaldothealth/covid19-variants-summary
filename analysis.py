import sys
import json
import argparse
import datetime
from pathlib import Path

import pandas as pd
from epiweeks import Week

# <- input
METADATA = "gisaid/metadata.tsv"
COUNTRIES = "data/countries.json"
VAX = "data/vaccinations.csv"

# -> output
COMPLETENESS = "completeness.csv"
GENDER = "gender.txt"
WEEKLY = "weekly.csv"
START_END_DATES = "start_end_dates.csv"

# other constants
VARIANTS = ["Omicron", "Delta"]


def _(s):
    print(">>>", s)


def w(s):
    "Get list from newline separated string"
    return [i.strip() for i in s.split("\n") if i.strip() != ""]


with open(COUNTRIES) as fp:
    COUNTRIES_MAP = json.load(fp)

COLS = w(
    """
Virus name
Type
Accession ID
Collection date
Location
Additional location information
Sequence length
Host
Patient age
Gender
Clade
Pango lineage
Pangolin version
Variant
AA Substitutions
Submission date
Is reference?
Is complete?
Is high coverage?
Is low coverage?
N-Content
GC-Content
"""
)

COMP_COLS = ["Country", "N"] + COLS

OUTCOLS = w(
    """
Country
Week
Omicron
Delta
Other_variants
Total
"""
)


def get_country(location):
    country_name = location.split("/")[1].strip()
    return COUNTRIES_MAP.get(country_name)


def get_location(location):
    location = location.split("/")
    country_name = location[1].strip()
    if len(location) == 3:
        return COUNTRIES_MAP.get(country_name), location[2].strip()
    else:
        return COUNTRIES_MAP.get(country_name), None


def is_float(n):
    try:
        float(n)
        return True
    except ValueError:
        return False


def to_epiweek(date):
    try:
        return str(Week.fromdate(datetime.date.fromisoformat(date)))
    except ValueError:
        return None


def check_collection_before_submission(df):
    _("Check that collection date is before submission")
    incorrect = df[df["Collection date"] > df["Submission date"]]
    if not incorrect.empty:
        print(incorrect)
    assert incorrect.empty
    return df


def completeness(df, country):
    country_df = df[df.Country == country]
    N = len(country_df)
    columns = set(COLS) - {"Patient age"}
    comp = {col: 100 * sum(~pd.isna(country_df[col])) / N for col in columns}
    comp.update(
        {
            "Patient age": 100 * country_df["Patient age"].map(is_float).sum() / N,
            "Country": country_df.iloc[0].Country,
            "N": N,
        }
    )
    return comp


def read_metadata(metadata):
    _(f"Read metadata < {metadata}")
    return pd.read_csv(metadata, dtype=str, delimiter="\t")


def read_vax(file):
    _(f"Read vaccination data < {file}")
    data = []
    iso_weeks = set()
    df = pd.read_csv(file).sort_values("date", ascending=False)
    df["week"] = df.date.map(to_epiweek)
    for v in df.itertuples():
        if (v.iso_code, v.week) not in iso_weeks:
            data.append(
                (
                    v.iso_code,
                    v.week,
                    v.people_vaccinated,
                    v.people_fully_vaccinated,
                    v.total_boosters,
                )
            )
            iso_weeks.add((v.iso_code, v.week))
    return pd.DataFrame(
        data=data,
        columns=["Country", "Week", "Vaccinated", "Fully_vaccinated", "Boosted"],
    )


def filter_human_hosts(df):
    _("Filter for human hosts")
    return df[df.Host == "Human"]


def parse_location(df):
    _(f"Parse location into country and region {len(df)}")
    # df[["Country", "Region"]] = pd.DataFrame(
    #     [get_location(item.Location) for item in df.itertuples()]
    # ^^ NOTE: This gives results that differ from the old country-only code!
    countries_regions = [get_location(item.Location) for item in df.itertuples()]
    df["Country"] = [x[0] for x in countries_regions]
    df["Region"] = [x[1] for x in countries_regions]
    return df


def filter_location(df, country=None, region=None):
    _(f"Filter for countries and regions under observation {len(df)}")
    if region is not None and country is None:
        raise ValueError(
            "If you specify region, country ISO3 code has to be specified as well"
        )
    if country is None:
        return df[df.Country.isin(COUNTRIES_MAP.values())]
    if region is None:
        return df[df.Country == country]
    else:
        return df[df.Region.isin(region.split(",")) & df.Country == country]


def calculate_epiweeks(df):
    _("Calculate epiweeks")
    df["Week"] = df["Collection date"].map(to_epiweek)
    return df


def calculate_completeness(df):
    _("Calculate completeness")
    return (
        pd.DataFrame([completeness(df, c) for c in COUNTRIES_MAP.values()])
        .reindex(columns=COMP_COLS)
        .set_index("Country")
    )


def gender_sets(df):
    _("Getting genders for each country")
    return [
        c + "\t" + "\t".join(sorted(set(df[df.Country == c].Gender.map(str))))
        for c in COUNTRIES_MAP.values()
    ]


def is_variant(variant):
    def func(v):
        return f"VOC {variant}" in v if isinstance(v, str) else False

    return func


def add_variant_columns(df):
    _("Adding variant columns")
    for v in VARIANTS:
        df[v] = df.Variant.map(is_variant(v))
    df["Other_variants"] = ~df[VARIANTS].any(axis=1)
    df["Total"] = df[VARIANTS + ["Other_variants"]].sum(axis=1)
    return df


def filter_valid_age(df, enabled=False):
    if not enabled:
        return df
    _("Filter data for entries with valid age")
    return df[df["Patient age"].map(is_float)]


def aggregate(df):
    _("Aggregate country data by epiweek")
    df = df[OUTCOLS]
    return df.groupby(["Country", "Week"]).agg("sum")


def merge_vax(df, vax, enabled=True):
    if not enabled:
        return df
    _("Merging with vaccination data from OWID")
    return df.merge(read_vax(vax), on=["Country", "Week"])


def get_start_dates(df, threshold_proportion=0.9):
    _("Getting start and end dates for variants")
    df = df.sort_values("Week")
    # start week is first time variant crossed 90%
    # end week is last time variant went below 90%
    for var in VARIANTS:
        df[f"{var}_prop"] = df[var] / df.Total
    data = []
    for country in sorted(set(df.Country)):
        cdf = df[df.Country == country]
        above_threshold = {
            var: cdf[cdf[f"{var}_prop"] >= threshold_proportion] for var in VARIANTS
        }
        below_threshold = {
            var: cdf[cdf[f"{var}_prop"] < threshold_proportion] for var in VARIANTS
        }

        start_dates = {
            var: str(above_threshold[var].iloc[0].Week)
            if not above_threshold[var].empty
            else None
            for var in VARIANTS
        }
        for var in sorted(VARIANTS):
            if start_dates[var] is None:
                end_date = None  # variant never reached threshold
            else:
                below_df = below_threshold[var]
                # end date should be after start date
                below_df = below_df[below_df.Week.astype(str) > start_dates[var]]
                end_date = str(below_df.iloc[0].Week) if not below_df.empty else None
            data.append(
                (country, var, start_dates[var], end_date, threshold_proportion)
            )
    return pd.DataFrame(
        data=data, columns=["Country", "Variant", "Start", "End", "Threshold"]
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        help="Input metadata file to use (default: gisaid/metadata.tsv)",
        default="gisaid/metadata.tsv",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output folder to use (default: output)",
        default="output",
    )
    parser.add_argument(
        "--valid-age",
        help="Filter to entries that have valid age",
        action="store_true",
    )
    parser.add_argument(
        "--region", help="Region to filter for, if specified, must have --country"
    )
    parser.add_argument("--country", help="ISO3 code of country to limit results to")
    args = parser.parse_args()
    OUTPUT = Path(args.output)
    if args.valid_age and not args.country:
        print("Filtering by age requires specifying country")
        sys.exit(1)
    if args.country:
        components = [args.country]
        if args.valid_age:
            components.append("valid-age")
        if args.region:
            components.extend(args.region.split(","))
        prefix = "_".join(components) + "_"
        # other outputs are not presented for country or region
        WEEKLY = prefix + WEEKLY

    if not OUTPUT.exists():
        OUTPUT.mkdir()

    df = (
        read_metadata(args.input)
        .pipe(filter_human_hosts)
        .pipe(parse_location)
        .pipe(filter_location, country=args.country, region=args.region)
        .pipe(filter_valid_age, enabled=args.valid_age)
        .pipe(check_collection_before_submission)
    )
    if args.country is None:
        calculate_completeness(df).to_csv(OUTPUT / COMPLETENESS)
        (OUTPUT / GENDER).write_text("\n".join(gender_sets(df)))
    # fmt: off
    weekly = (
        df
        .pipe(calculate_epiweeks)
        .pipe(add_variant_columns)
        .pipe(aggregate)
        .pipe(merge_vax, vax=VAX, enabled=(not args.country))
    )
    # fmt: on
    weekly.to_csv(OUTPUT / WEEKLY, index=bool(args.country))
    if args.country is None:
        get_start_dates(weekly).to_csv(OUTPUT / START_END_DATES, index=False)
