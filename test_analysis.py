import pandas as pd
import pytest

import analysis


@pytest.mark.parametrize(
    "source,expected",
    [
        ("South America / Argentina / Buenos Aires", ("ARG", "Buenos Aires")),
        ("South America / Argentina", ("ARG", None)),
    ],
)
def test_get_location(source, expected):
    assert analysis.get_location(source) == expected


def test_calculate_epiweeks():
    dates = {"Collection date": ["2022-01-01", "2022-01-20", "2021-12-12"]}
    weeks = {"Week": ["202152", "202203", "202150"]}
    df = pd.DataFrame({**dates})
    expected_df = pd.DataFrame({**dates, **weeks})
    assert analysis.calculate_epiweeks(df).equals(expected_df)
