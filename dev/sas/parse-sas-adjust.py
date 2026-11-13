# Parse the Tukey-Kramer differences of pm-reference-adjust.sas (results page as text,
# pdftotext -layout) into a CSV: python3 -I dev/sas/parse-sas-adjust.py <text> <csv>
import sys, csv, re
txt, out = sys.argv[1], sys.argv[2]
lines = open(txt, encoding="utf-8", errors="replace").read().splitlines()
run = 0
rows = []
for ln in lines:
    if "The Mixed Procedure" in ln:
        run += 1
    tok = ln.split()
    if not tok or tok[0] not in ("Variety", "Variety*Nitro"):
        continue
    if "Tukey-" not in tok:
        continue            # LS-means rows
    i = tok.index("Tukey-")
    if tok[0] == "Variety":
        lev = tok[1:3]; lev = [lev[0], "", lev[1], ""]; num = tok[3:i]
    else:
        lev = tok[1:5]; num = tok[5:i]
    est, se, df, t, p = num
    adjp, alpha, lo, up, alo, aup = tok[i + 1:i + 7]
    rows.append(dict(adjdfe="row" if run == 1 else "source", effect=tok[0], variety=lev[0], nitro=lev[1],
                     variety2=lev[2], nitro2=lev[3], estimate=est, se=se, df=df, t=t, p=p, adj_p=adjp,
                     lower=lo, upper=up, adj_lower=alo, adj_upper=aup))
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
print(len(rows), "rows;", sum(r["adjdfe"] == "row" for r in rows), "row,", sum(r["adjdfe"] == "source" for r in rows), "source")
