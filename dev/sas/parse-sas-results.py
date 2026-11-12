"""Parse the text of a SAS Studio results page of pm-reference.sas (printed to
PDF, then `pdftotext -layout`) into one CSV with a row per table row.

Usage: python3 -I dev/sas/parse-sas-results.py results.txt out.csv

Columns: run (1-16), data, ddfm, table, effect, level1, level2, level3, level4
(the class levels; for differences the first and the second mean), and the
numbers of the row (est, se, df, t, p, lower, upper for LSMEANS and
differences; numdf, dendf, F, p for Type 3 tests; est for covariance
parameters; est, se, df, t, p for fixed-effect solutions). Values keep the
precision SAS printed.
"""
import csv
import re
import sys

src, out = sys.argv[1], sys.argv[2]
lines = open(src, encoding="utf-8").read().split("\n")

TITLES = ["Covariance Parameter Estimates", "Fit Statistics", "Solution for Fixed Effects",
          "Type 3 Tests of Fixed Effects", "Least Squares Means", "Differences of Least Squares Means",
          "Iteration History", "Number of Observations", "Dimensions", "Class Level Information",
          "Model Information", "Convergence criteria met"]
EFFECTS = {"Intercept", "Variety", "Nitro", "Variety*Nitro", "Group", "Time", "Group*Time"}
DDFM = {"Containment": "CONTAIN", "Between-Within": "BETWITHIN", "Satterthwaite": "SATTERTHWAITE",
        "Kenward-Roger": "KR"}


def num(x):
    if x == ".":
        return ""
    if x.startswith("<"):
        return x            # e.g. <.0001
    return x


def is_num(x):
    return bool(re.fullmatch(r"-?\d*\.?\d+(E[-+]?\d+)?|<\.\d+|\.", x))


rows = []
run = 0
data = ddfm = table = None
for line in lines:
    s = line.strip()
    if not s or s.startswith("https://") or "Results: Program 1" in s:
        continue
    m = re.match(r"Data Set\s+WORK\.(\S+)", s)
    if m:
        run += 1
        data = m.group(1)
        table = None
        continue
    m = re.match(r"Number of Observations Used\s+(\d+)", s)
    if m:
        rows.append(dict(run=run, data=data, ddfm=ddfm, table="nobs", effect="", est=m.group(1)))
        continue
    m = re.match(r"Degrees of Freedom Method\s+(.+)$", s)
    if m:
        ddfm = DDFM[m.group(1).strip()]
        continue
    hit = [t for t in TITLES if s.startswith(t)]
    if hit:
        table = hit[0]
        continue
    tok = s.split()
    if table == "Covariance Parameter Estimates":
        if tok[0] in ("Cov",):
            continue
        if is_num(tok[-1]) and not is_num(tok[0]):
            rows.append(dict(run=run, data=data, ddfm=ddfm, table="covparms", effect=tok[0],
                             level1=tok[1] if len(tok) == 3 else "", est=tok[-1]))
        continue
    if tok[0] not in EFFECTS:
        continue
    if table == "Type 3 Tests of Fixed Effects" and len(tok) == 5:
        rows.append(dict(run=run, data=data, ddfm=ddfm, table="tests3", effect=tok[0], numdf=tok[1],
                         dendf=tok[2], F=tok[3], p=num(tok[4])))
    elif table == "Solution for Fixed Effects":
        # Effect levels... Estimate SE DF t p  (aliased rows: 0 . . . .)
        vals, levs = tok[-5:], tok[1:-5]
        rows.append(dict(run=run, data=data, ddfm=ddfm, table="solution", effect=tok[0],
                         **{f"level{i + 1}": v for i, v in enumerate(levs)},
                         est=vals[0], se=num(vals[1]), df=num(vals[2]), t=num(vals[3]), p=num(vals[4])))
    elif table in ("Least Squares Means", "Differences of Least Squares Means"):
        # Effect levels... Estimate SE DF t p Alpha Lower Upper
        vals, levs = tok[-8:], tok[1:-8]
        if not all(is_num(v) for v in vals):
            continue
        rows.append(dict(run=run, data=data, ddfm=ddfm,
                         table="lsmeans" if table == "Least Squares Means" else "diffs", effect=tok[0],
                         **{f"level{i + 1}": v for i, v in enumerate(levs)},
                         est=vals[0], se=vals[1], df=vals[2], t=vals[3], p=num(vals[4]), lower=vals[6],
                         upper=vals[7]))

fields = ["run", "data", "ddfm", "table", "effect", "level1", "level2", "level3", "level4", "est", "se", "df",
          "t", "p", "lower", "upper", "numdf", "dendf", "F"]
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    for r in rows:
        w.writerow(r)
print(f"{run} runs, {len(rows)} rows")
