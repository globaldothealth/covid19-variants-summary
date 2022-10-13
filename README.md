# COVID-19 variants summary

This code summarises data from GISAID by epiweek. Outputs and analyses using these data are published in Gonçalves et al. eLife 2022;11:e80556. DOI: https://doi.org/10.7554/eLife.80556 

[![tests](https://github.com/globaldothealth/covid19-variants-summary/actions/workflows/tests.yml/badge.svg)](https://github.com/globaldothealth/covid19-variants-summary/actions/workflows/tests.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7195830.svg)](https://doi.org/10.5281/zenodo.7195830)


## Requirements

Python 3.6 or later.

Python modules: [pandas](https://pandas.pydata.org),
[epiweeks](https://epiweeks.readthedocs.io)

To install the modules in a virtual environment:

    python3 -m venv venv
    source venv/bin/activate
    pip install -r requirements.txt
    python3 analysis.py  # run analysis


## Data

Place GISAID metadata at *gisaid/metadata.tsv*

## Output

All output files are in the *output* folder.

* **weekly.csv**: Main output file with weekly totals of VOC Omicron and
  Delta as well as other variants by country. We use the epiweeks
  package to calculate CDC epiweeks (starting on Sunday).

* **completeness.csv**: Completeness as percentage of fields which are
  not null, other than the Patient age field which is the percentage of
  fields which are of type float.

* **gender.txt**: The gender field in GISAID is particularly error-prone
  with many instances of ages and mis-spellings. This lists all the
  unique values in the Gender field by country.
