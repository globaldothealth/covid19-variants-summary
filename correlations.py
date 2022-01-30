import json
import contextlib
from pathlib import Path
import pandas as pd

from analysis import VARIANTS, COUNTRIES_MAP

CORRELATIONS = "output/correlations.json"


def correlations(df_index, df_other):
    df = df_index.merge(df_other, on=["Country", "Week"])
    return {
        variant: df[variant + "_x"].corr(df[variant + "_y"])
        for variant in VARIANTS + ["Total"]
    }


def subset_correlations(country, output="output"):
    files = list(Path.cwd().glob(f"{output}/{country}/{country}*.csv"))
    index_file_name = f"{country}_weekly.csv"
    if not any(f.name == index_file_name for f in files):
        raise ValueError(f"Generate index file {output}/{country}/{country}_weekly.csv")
    df_index = pd.read_csv(f"{output}/{country}/{index_file_name}")
    total = df_index.Total.sum()
    retval = {}
    for file in [
        f for f in files if f.name != index_file_name and "group" not in f.name
    ]:
        df = pd.read_csv(file)
        corr = correlations(df_index, df)
        corr["ratio"] = df.Total.sum() / total
        retval[file.name] = corr
    return retval


def week_to_float(weekstr: str) -> float:
    w = str(weekstr)
    assert len(w) == 6
    year, week = w[:4], w[4:]
    return int(year) + (int(week) - 1) / 52


_correlations = {}

for country in COUNTRIES_MAP.values():
    with contextlib.suppress(ValueError):
        _correlations[country] = subset_correlations(country)
with open(CORRELATIONS, "w") as fp:
    json.dump(_correlations, fp, indent=2, sort_keys=True)
