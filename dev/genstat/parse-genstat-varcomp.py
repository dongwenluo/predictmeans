## Variance components of dev/genstat/pm-reference-3.gen (output as text:
## genstat-results-2026-11-14-run3.txt) to CSV: case, term, est, se (as
## printed). The REML analyses are numbered in the order of the program.
## Usage: python3 parse-genstat-varcomp.py results.txt out.csv
import sys, csv, re
cases = ["DYESTUFF", "DYESTUFF2", "PENICILLIN", "ORTHO_SLOPE_CENTRED", "ORTHO_AR1"]
num = re.compile(r"^-?[0-9.]+$")
rows, k, sec = [], -1, None
for line in open(sys.argv[1]):
    s = line.strip()
    c = s.split()
    if s.startswith("REML variance components analysis"):
        k += 1
        sec = None
    elif s.startswith("Estimated variance components"):
        sec = "vc"
    elif s.startswith("Residual variance model"):
        sec = "res"
    elif sec == "vc" and len(c) == 3 and num.match(c[1]) and num.match(c[2]):
        rows.append([cases[k], c[0], c[1], c[2]])
    elif sec == "res" and len(c) >= 3 and c[-3] in ("Sigma2", "phi_1") and num.match(c[-2]):
        rows.append([cases[k], c[-3], c[-2], c[-1]])
with open(sys.argv[2], "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["case", "term", "est", "se"])
    w.writerows(rows)
print(len(rows), "rows")
