import random
import argparse
import datetime

TODAY = datetime.date.today()


def virus_name():
    while True:
        yield f"hCoV-19/Argentina/FAKE-{random.choice(range(200000, 300000))}/2021"


def accession_id():
    n = 199999
    while True:
        n += 1
        yield f"EPI_FAKE_{n}"


def from_range(start, end):
    while True:
        yield random.choice(range(start, end))


def pick(population, weights=None):
    while True:
        yield random.choices(population, weights)[0]


def from_date(start: datetime.date, end: datetime.date = None):
    days = ((end or TODAY) - start).days
    while True:
        yield str(start + datetime.timedelta(days=random.choice(range(days))))


def const(x):
    while True:
        yield x


def true(probability):
    "Yield True with a certain probability"

    while True:
        yield str(random.random() > probability)


COLS = [
    ("Virus name", virus_name()),
    ("Type", "betacoronavirus"),
    ("Accession ID", accession_id()),
    ("Collection date", from_date(datetime.date(2021, 6, 1))),
    ("Location", "South America / Argentina / Buenos Aires"),
    ("Additional location information", ""),
    ("Sequence length", from_range(29000, 30000)),
    ("Host", pick(["Human", "Mouse"], [90, 10])),
    ("Patient age", from_range(0, 100)),
    ("Gender", pick(["Male", "Female", "Other"], [0.3, 0.3, 0.2])),
    ("Clade", "QQ"),
    ("Pango lineage", "ZZ.9"),
    ("Pangolin version", "2021-12-06"),
    ("Variant", pick(["VOC Omicron", "VOC Delta"], [0.4, 0.6])),
    ("AA Substitutions", "AA"),
    ("Submission date", str(TODAY)),
    ("Is reference?", ""),
    ("Is complete?", "True"),
    ("Is high coverage?", true(0.75)),
    ("Is low coverage?", true(0.2)),
    ("N-Content", "10"),
    ("GC-Content", "55"),
]


def generate(N):
    "Generate N entries"
    print("\t".join(col for col, _ in COLS))
    for _ in range(N):
        print("\t".join(gen if isinstance(gen, str) else str(next(gen)) for _, gen in COLS))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-n",
        "--number",
        help="How many entries to generate (default: 100)",
        default=100,
        type=int,
    )
    parser.add_argument(
        "-s", "--seed", help="Random number seed to use (default: 10)", default=10
    )
    args = parser.parse_args()
    random.seed(args.seed)
    generate(args.number)
