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


@pytest.mark.parametrize(
    "source,expected",
    [
        ({"Patient age lower": "20", "Patient age upper": "25"}, "20 – 30"),
        ({"Patient age lower": "3", "Patient age upper": "4"}, "1 – 5"),
        ({"Patient age lower": "20", "Patient age upper": "35"}, ""),
    ],
)
def test_get_age_group_row(source, expected):
    assert analysis.get_age_group_row(source) == expected


@pytest.mark.parametrize(
    "source,expected",
    [
        ({"Patient age": "20", "Gender": "Male"}, ("20", "20", "Male")),
        ({"Patient age": 20, "Gender": "Male"}, (20, 20, "Male")),
        ({"Patient age": "Male", "Gender": "20"}, ("20", "20", "Male")),
        ({"Patient age": "20 - 29", "Gender": "F"}, ("20.0", "29.0", "F")),
        ({"Patient age": "29 - 20", "Gender": "F"}, ("", "", "F")),
        ({"Patient age": "30 to 35", "Gender": "F"}, ("30.0", "35.0", "F")),
        ({"Patient age": "F", "Gender": "30 to 35"}, ("30.0", "35.0", "F")),
        ({"Patient age": "30 - 50", "Gender": "M"}, ("30.0", "50.0", "M")),
    ],
)
def test_fix_age_gender_row(source, expected):
    assert analysis.fix_age_gender_row(source) == expected


def test_calculate_epiweeks():
    dates = {"Collection date": ["2022-01-01", "2022-01-20", "2021-12-12"]}
    weeks = {"Week": ["202152", "202203", "202150"]}
    df = pd.DataFrame({**dates})
    expected_df = pd.DataFrame({**dates, **weeks})
    assert analysis.calculate_epiweeks(df).equals(expected_df)
