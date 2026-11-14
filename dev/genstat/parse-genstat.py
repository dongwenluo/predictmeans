## Genstat REML / VPREDICT / ANOVA output (text of the RTF) -> CSV rows:
## analysis, table, term, key, value (as printed). Usage: parse.py run1.txt run2.txt out.csv
import sys, csv, re
rows = []
num = re.compile(r"^-?[\d.]+$")
def add(an, table, term, key, value):
    rows.append(dict(analysis=an, table=table, term=term, key=key, value=value.strip()))
def cells(line):
    return [x.strip() for x in line.split("\t") if x.strip() != ""]
for run, path in enumerate(sys.argv[1:-1], start=1):
    L = [l.rstrip("\n") for l in open(path)]
    an = sec = sub = term = None
    header = None; same = []; i = 0
    while i < len(L):
        l = L[i]; s = l.strip(); c = cells(l)
        if s.startswith("REML variance components analysis"): sec, sub = "reml", None
        elif s.startswith("Predictions from REML analysis"): sec, sub, term = "vpredict", None, None
        elif s.startswith("Analysis of variance"): sec, sub = "anova", None
        elif s.startswith("Tables of means"): sec, sub, term = "anova_means", None, None
        elif s.startswith("Response variate:") or (s.startswith("Variate:") and sec.startswith("anova")):
            an = f"run{run}:{s.split()[-1]}" + (":anova" if sec.startswith("anova") else "")
        elif sec == "reml":
            if c and c[0] in ("Block", "Block.Variety") and len(c) == 3 and num.match(c[1]):
                add(an, "vc", c[0], "component", c[1])
                add(an, "vc", c[0], "se", c[2])
            elif len(c) >= 4 and c[-3] in ("Sigma2", "phi_1"):
                add(an, "vc", c[-3], "estimate", c[-2])
                add(an, "vc", c[-3], "se", c[-1])
            elif s.startswith("Sequentially adding terms"): sub = "seq"
            elif s.startswith("Dropping individual terms"): sub = "drop"
            elif sub in ("seq", "drop") and len(c) == 6 and num.match(c[1]):
                for k, v in zip(("wald", "ndf", "F", "ddf", "p"), c[1:]): add(an, "wald_" + sub, c[0], k, v)
            elif s.startswith("Table of predicted means for"):
                term = s.split("for", 1)[1].strip(); sub = "means"; header = None
            elif sub == "means" and header is None and len(c) >= 3:
                header = c[1:] if "." not in term else c[1:]
                if "." in term:                       # two-way: "Nitro <levels>", then "Variety", then rows
                    sub = "means2"
            elif sub == "means" and header is not None and c and num.match(c[0]):
                for lev, v in zip(header, c): add(an, "mean", term, lev, v)
                sub = "sed"
            elif sub == "means2" and len(c) == len(header) + 1:
                for lev, v in zip(header, c[1:]): add(an, "mean", term, f"{c[0]}:{lev}", v)
            elif s.startswith("Standard error of differences:"):
                add(an, "sed", term, "single", s.split(":")[1])
            elif s.startswith("Standard error of differences for same level"):
                same = cells(L[i + 1]) if cells(L[i + 1]) else cells(L[i + 2]); sub = "same"
            elif c and c[0] in ("Average:", "Maximum:", "Minimum:"):
                k = c[0][:-1].lower()
                if sub == "same" and len(c) == 3:
                    add(an, "sed_same", term, f"{k}:{same[0]}", c[1]); add(an, "sed_same", term, f"{k}:{same[1]}", c[2])
                elif len(c) == 2:
                    add(an, "sed", term, k, c[1])
        elif sec == "vpredict":
            if len(c) >= 3 and c[2] == "Classifies predictions": term = c[0]
            elif s == "Predictions": sub, header = "vp_pred", None
            elif s == "Standard errors": sub, header = "vp_se", None
            elif s == "Standard error of differences": sub, header = "vp_sed", []
            elif sub in ("vp_pred", "vp_se") and c:
                if header is None and c[0] == term: header = c[1:]
                elif header is not None and num.match(c[0]):
                    for lev, v in zip(header, c): add(an, sub, term, lev, v)
                    sub = None
            elif sub == "vp_sed" and len(c) >= 3 and c[0].startswith(term + " "):
                lev = c[0][len(term) + 1:]; k = int(c[1]); header.append(lev)
                for j, v in enumerate(c[2:2 + k]):
                    if v != "*": add(an, "vp_sed", term, f"{header[j]}|{lev}", v)
        elif sec == "anova":
            if c and c[0] in ("Variety", "Residual", "Nitro", "Variety.Nitro", "Total") and len(c) >= 3 and re.match(r"^\d+$", c[1]):
                add(an, "aov_table", c[0], "df|...", "|".join(c[1:]))
        elif sec == "anova_means":
            if s.startswith("Grand mean"): add(an, "aov_mean", "Grand", "mean", s.split()[-1])
            elif c and c[0] in ("Variety", "Nitro") and len(c) > 2 and c[1] != "Nitro":
                header = c[1:]; vals = cells(L[i + 1])
                for lev, v in zip(header, vals): add(an, "aov_mean", c[0], lev, v)
                i += 2; continue
            elif len(c) > 3 and c[0] == "Variety" and c[1] == "Nitro":
                header = c[2:]; j = i + 1
                while j < len(L) and cells(L[j]):
                    v = cells(L[j])
                    for lev, x in zip(header, v[1:]): add(an, "aov_mean", "Variety.Nitro", f"{v[0]}:{lev}", x)
                    j += 1
                i = j; continue
            elif c and c[0] in ("rep.", "s.e.d.", "d.f."): add(an, "aov_sed", c[0], "values", "|".join(c[1:]))
            elif len(c) == 2 and c[0] == "Variety" and num.match(c[1]): add(an, "aov_sed", "same Variety", "s.e.d.", c[1])
            elif s.startswith("Unit"): sub = "mv"
            elif sub == "mv" and len(c) == 2 and re.match(r"^\d+$", c[0]):
                add(an, "missing_value", "unit " + c[0], "estimate", c[1]); sub = None
        i += 1
with open(sys.argv[-1], "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["analysis", "table", "term", "key", "value"]); w.writeheader(); w.writerows(rows)
print(len(rows), "rows")
