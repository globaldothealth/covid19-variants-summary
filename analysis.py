import json
import argparse
import datetime
from pathlib import Path

import pandas as pd
from epiweeks import Week

# <- input
METADATA = "gisaid/metadata.tsv"
COUNTRIES = "data/countries.json"

# -> output
COMPLETENESS = "completeness.csv"
GENDER = "gender.txt"
WEEKLY = "weekly.csv"

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
    print(country, region)
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


def aggregate(df):
    _("Aggregate country data by epiweek")
    df = df[OUTCOLS]
    return df.groupby(["Country", "Week"]).agg("sum")


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
        default="output"
    )
    parser.add_argument(
        "--region", help="Region to filter for, if specified, must have --country"
    )
    parser.add_argument("--country", help="ISO3 code of country to limit results to")
    args = parser.parse_args()
    OUTPUT = Path(args.output)
    if not OUTPUT.exists():
        OUTPUT.mkdir()

    df = (
        read_metadata(args.input)
        .pipe(filter_human_hosts)
        .pipe(parse_location)
        .pipe(filter_location, country=args.country, region=args.region)
        .pipe(check_collection_before_submission)
    )
    calculate_completeness(df).to_csv(OUTPUT / COMPLETENESS)
    (OUTPUT / GENDER).write_text("\n".join(gender_sets(df)))
    # fmt: off
    (
        df
        .pipe(calculate_epiweeks)
        .pipe(add_variant_columns)
        .pipe(aggregate)
    ).to_csv(OUTPUT / WEEKLY)
    # fmt: on
