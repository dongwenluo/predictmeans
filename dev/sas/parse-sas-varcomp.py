## Covariance parameter tables of dev/sas/pm-reference-varcomp.sas (results
## page as text: sas-results-2026-11-14-varcomp.txt) to CSV: case, parm,
## subject, est, se, lower, upper (as printed; "." for none).
## Usage: python3 parse-sas-varcomp.py results.txt out.csv
import sys, csv, re
cases = [("Oats_MissWP", "OATS_MISSWP"), ("Oats_Miss ", "OATS_MISS"), ("Oats (random", "OATS"),
         ("intercept and slope", "ORTHO_SLOPE"), ("intercept and AR(1)", "ORTHO_AR1"),
         ("unstructured repeated", "ORTHO_UN"), ("AR(1) repeated measures", "AR1"),
         ("Dyestuff2", "DYESTUFF2"), ("Dyestuff (random", "DYESTUFF"), ("Penicillin", "PENICILLIN")]
num = re.compile(r"^-?[0-9.]+(E-?[0-9]+)?$|^\.$")
rows, case, inblock = [], None, False
for line in open(sys.argv[1]):
    s = line.strip()
    if "varcomp:" in s:
        case = next(k for t, k in cases if t in s + " ")
        continue
    if s.startswith("Covariance Parameter Estimates"):
        inblock = True
        continue
    if s.startswith("Asymptotic Covariance") or s.startswith("Fit Statistics"):
        inblock = False
        continue
    if not inblock:
        continue
    t = s.split()
    if len(t) == 9 and all(num.match(x) or x == "<.0001" for x in t[2:]):
        parm, subj, rest = t[0], t[1], t[2:]
    elif len(t) == 8 and all(num.match(x) or x == "<.0001" for x in t[1:]):
        parm, subj, rest = t[0], "", t[1:]
    else:
        continue
    rows.append([case, parm, subj, rest[0], rest[1], rest[5], rest[6]])
with open(sys.argv[2], "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["case", "parm", "subject", "est", "se", "lower", "upper"])
    w.writerows(rows)
print(len(rows), "rows")
