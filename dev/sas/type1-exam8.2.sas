/* Type 1 (sequential) F tests in a mixed model: which hypotheses SAS uses.
   For predictmeans::anova(type = "I") (R/anova.R, pm_type1_L()).

   Data: Example 8.2 of Williams, Matheson and Harwood (2002), Eucalyptus
   clone trial (eda4treeR::DataExam8.2, public): 5 replicates x 10 rows x 6
   columns, 56 test clones and 4 standards (standards have clone "0").
   Columns: repl row column contcompf (1 = test clone, 0 = standard)
   std (standard name, "0" for test clones) clone cv (= contcompf, numeric)
   dbh.

   Model: dbh ~ repl + column + contcompf + contcompf:std, random
   repl:row, repl:column and clone effects for test clones only (cv | clone,
   = Genstat contcompv.clone with contcompv 0/1), REML.

   Question: Genstat's "sequentially adding terms" Wald statistics come from
   the GLS fit (X'V^-1 X); lmerTest builds Type I hypotheses from X'X. They
   agree for the first and the last term, and differ for column and contcompf:

                      Genstat / GLS        X'X (lmerTest)
     term             Wald    F            F (unadjusted)
     repl             8.20    2.049        2.049
     column           7.99    1.599        1.487
     contcompf       72.66   72.66        74.49
     contcompf*std   51.75   17.25        17.25

   Genstat d.d.f. (KR, average information): 26.6 19.7 254.0 208.2.
   R (lme4 REML) variance components: repl:row 0.0802, repl:column 0.0529,
   clone 0.4950, residual 0.3992.

   What to send back: the Tests1 and Tests3 CSV files of both runs (or the
   results page). Run 1 is DDFM=KR (observed information), run 2 DDFM=SATTERTHWAITE.
   If SAS's Type 1 F for column is 1.60 (KR may scale it slightly), SAS uses
   the GLS hypotheses; if 1.49, X'X. */

data ex82;
  input repl row column contcompf std $ clone $ cv dbh;
  datalines;
1 1 1 1 0 UP188 1 10.33
1 1 2 0 UG323 0 0  9.27
1 1 3 1 0 UP266 1  8.86
1 1 4 1 0 UP227 1  9.14
1 1 5 1 0 UP254 1 10.11
1 1 6 1 0 UP110 1  7.00
1 2 1 1 0 UP219 1  9.30
1 2 2 1 0 UP297 1  8.28
1 2 3 1 0 UP278 1  8.38
1 2 4 1 0 UP287 1  7.91
1 2 5 1 0 UP115 1  7.26
1 2 6 1 0 UP150 1  7.49
1 3 1 1 0 UP316 1  7.06
1 3 2 1 0 UP218 1  8.07
1 3 3 1 0 UP147 1  8.40
1 3 4 1 0 UP276 1  9.31
1 3 5 1 0 UP111 1  7.35
1 3 6 1 0 UP240 1  8.30
1 4 1 1 0 UP220 1  8.96
1 4 2 1 0 UP251 1  8.23
1 4 3 1 0 UP288 1  8.97
1 4 4 1 0 UP156 1  6.60
1 4 5 1 0 UP192 1  8.47
1 4 6 1 0 UP144 1  7.30
1 5 1 1 0 UP171 1 10.27
1 5 2 1 0 UP159 1  7.87
1 5 3 1 0 UP317 1  7.09
1 5 4 1 0 UP124 1  8.25
1 5 5 1 0 UP274 1 10.83
1 5 6 1 0 UP301 1  9.19
1 6 1 1 0 UP151 1  8.80
1 6 2 1 0 UP138 1  8.97
1 6 3 1 0 UP291 1  8.96
1 6 4 1 0 UP296 1  8.34
1 6 5 0 U6 0 0  6.89
1 6 6 1 0 UP122 1 10.57
1 7 1 1 0 UP300 1  9.85
1 7 2 1 0 UP283 1  8.52
1 7 3 1 0 UP149 1  8.89
1 7 4 1 0 UP165 1  8.27
1 7 5 1 0 UP233 1  6.96
1 7 6 1 0 UP169 1  7.83
1 8 1 1 0 UP123 1  8.92
1 8 2 1 0 UP153 1  9.78
1 8 3 1 0 UP282 1  8.64
1 8 4 1 0 UP187 1  7.41
1 8 5 1 0 UP140 1  8.21
1 8 6 1 0 UP127 1  7.24
1 9 1 1 0 UP135 1  9.06
1 9 2 1 0 UP190 1  8.65
1 9 3 1 0 UP262 1  9.00
1 9 4 1 0 UP223 1  9.74
1 9 5 1 0 UP225 1  9.07
1 9 6 1 0 UP238 1  8.88
1 10 1 1 0 UP290 1  7.97
1 10 2 0 PN14 0 0  6.52
1 10 3 1 0 UP154 1  8.14
1 10 4 1 0 UP103 1  8.25
1 10 5 0 SSOseed 0 0  6.28
1 10 6 1 0 UP252 1  9.05
2 1 1 1 0 UP297 1  9.22
2 1 2 1 0 UP288 1  9.40
2 1 3 1 0 UP140 1  9.62
2 1 4 1 0 UP138 1  9.20
2 1 5 1 0 UP254 1  9.50
2 1 6 1 0 UP290 1  8.35
2 2 1 1 0 UP103 1  8.52
2 2 2 1 0 UP123 1  8.22
2 2 3 0 PN14 0 0  7.11
2 2 4 1 0 UP144 1  8.24
2 2 5 1 0 UP291 1  8.31
2 2 6 1 0 UP276 1  9.91
2 3 1 1 0 UP316 1  8.07
2 3 2 1 0 UP251 1  9.61
2 3 3 1 0 UP124 1  8.21
2 3 4 0 SSOseed 0 0  5.74
2 3 5 1 0 UP300 1  9.04
2 3 6 1 0 UP190 1  7.82
2 4 1 1 0 UP266 1  9.47
2 4 2 1 0 UP240 1  8.86
2 4 3 1 0 UP169 1  8.38
2 4 4 1 0 UP135 1  9.75
2 4 5 1 0 UP159 1  7.61
2 4 6 1 0 UP220 1  9.15
2 5 1 1 0 UP122 1  8.90
2 5 2 1 0 UP150 1  6.61
2 5 3 1 0 UP192 1  8.79
2 5 4 1 0 UP188 1  9.35
2 5 5 1 0 UP301 1  9.43
2 5 6 1 0 UP165 1  8.83
2 6 1 1 0 UP262 1  9.54
2 6 2 1 0 UP115 1  8.94
2 6 3 0 UG323 0 0  9.38
2 6 4 0 U6 0 0  5.70
2 6 5 1 0 UP283 1  9.04
2 6 6 1 0 UP317 1  6.70
2 7 1 1 0 UP153 1 10.71
2 7 2 1 0 UP219 1  8.00
2 7 3 1 0 UP171 1  8.91
2 7 4 1 0 UP238 1  9.30
2 7 5 1 0 UP252 1  9.10
2 7 6 1 0 UP111 1  8.18
2 8 1 1 0 UP156 1  7.54
2 8 2 1 0 UP282 1 10.09
2 8 3 1 0 UP223 1  9.62
2 8 4 1 0 UP218 1  9.90
2 8 5 1 0 UP274 1  9.99
2 8 6 1 0 UP154 1  8.37
2 9 1 1 0 UP233 1  8.43
2 9 2 1 0 UP227 1  9.12
2 9 3 1 0 UP278 1  9.18
2 9 4 1 0 UP147 1  9.74
2 9 5 1 0 UP151 1 10.14
2 9 6 1 0 UP127 1  7.70
2 10 1 1 0 UP225 1  7.94
2 10 2 1 0 UP296 1  9.73
2 10 3 1 0 UP149 1 10.04
2 10 4 1 0 UP287 1 10.16
2 10 5 1 0 UP187 1  8.00
2 10 6 1 0 UP110 1  6.85
3 1 1 1 0 UP190 1  8.47
3 1 2 1 0 UP122 1 10.18
3 1 3 1 0 UP115 1  9.01
3 1 4 1 0 UP187 1  7.21
3 1 5 1 0 UP147 1  9.07
3 1 6 1 0 UP171 1  8.08
3 2 1 1 0 UP283 1  9.19
3 2 2 1 0 UP240 1  8.06
3 2 3 0 SSOseed 0 0  6.20
3 2 4 1 0 UP103 1  8.74
3 2 5 1 0 UP278 1  9.51
3 2 6 1 0 UP220 1  9.22
3 3 1 1 0 UP266 1  9.62
3 3 2 1 0 UP138 1  9.41
3 3 3 1 0 UP223 1  9.11
3 3 4 1 0 UP219 1  8.16
3 3 5 1 0 UP127 1  7.52
3 3 6 1 0 UP274 1  9.67
3 4 1 1 0 UP140 1  8.02
3 4 2 1 0 UP300 1 11.73
3 4 3 1 0 UP297 1  9.74
3 4 4 1 0 UP317 1  9.21
3 4 5 1 0 UP296 1 10.22
3 4 6 1 0 UP144 1  8.02
3 5 1 1 0 UP123 1  9.72
3 5 2 1 0 UP124 1  8.42
3 5 3 1 0 UP110 1  6.99
3 5 4 1 0 UP218 1  8.40
3 5 5 1 0 UP254 1  8.85
3 5 6 1 0 UP262 1  9.68
3 6 1 1 0 UP276 1 10.54
3 6 2 1 0 UP233 1  9.89
3 6 3 1 0 UP282 1 10.78
3 6 4 1 0 UP252 1  9.46
3 6 5 1 0 UP135 1  8.70
3 6 6 0 UG323 0 0  8.82
3 7 1 1 0 UP288 1 10.69
3 7 2 0 PN14 0 0  7.86
3 7 3 1 0 UP151 1  8.45
3 7 4 1 0 UP188 1  8.56
3 7 5 1 0 UP301 1  8.63
3 7 6 1 0 UP225 1  8.53
3 8 1 1 0 UP316 1  8.99
3 8 2 1 0 UP154 1  8.87
3 8 3 1 0 UP169 1  8.02
3 8 4 1 0 UP192 1  7.88
3 8 5 1 0 UP153 1  9.60
3 8 6 1 0 UP227 1  8.46
3 9 1 1 0 UP238 1  9.53
3 9 2 1 0 UP287 1  9.76
3 9 3 1 0 UP165 1  9.64
3 9 4 1 0 UP251 1 10.32
3 9 5 1 0 UP111 1  7.30
3 9 6 1 0 UP291 1  9.43
3 10 1 0 U6 0 0  7.18
3 10 2 1 0 UP290 1  8.62
3 10 3 1 0 UP159 1  7.34
3 10 4 1 0 UP149 1  8.68
3 10 5 1 0 UP156 1  6.38
3 10 6 1 0 UP150 1  7.21
4 1 1 1 0 UP276 1 10.08
4 1 2 1 0 UP278 1  9.35
4 1 3 1 0 UP296 1  9.75
4 1 4 1 0 UP159 1  7.40
4 1 5 1 0 UP188 1  8.68
4 1 6 1 0 UP190 1  8.83
4 2 1 1 0 UP123 1 10.22
4 2 2 1 0 UP122 1  9.03
4 2 3 1 0 UP147 1  9.80
4 2 4 1 0 UP144 1  8.69
4 2 5 1 0 UP149 1  9.47
4 2 6 0 UG323 0 0  9.40
4 3 1 0 SSOseed 0 0  6.38
4 3 2 1 0 UP227 1 10.19
4 3 3 1 0 UP192 1  9.48
4 3 4 1 0 UP135 1  9.55
4 3 5 1 0 UP291 1 10.11
4 3 6 1 0 UP219 1  9.87
4 4 1 1 0 UP169 1  9.29
4 4 2 0 U6 0 0  7.52
4 4 3 0 PN14 0 0  8.97
4 4 4 1 0 UP251 1  9.60
4 4 5 1 0 UP127 1  7.06
4 4 6 1 0 UP262 1  9.74
4 5 1 1 0 UP187 1  7.56
4 5 2 1 0 UP288 1 10.87
4 5 3 1 0 UP274 1 11.65
4 5 4 1 0 UP165 1  8.46
4 5 5 1 0 UP297 1  9.31
4 5 6 1 0 UP240 1  7.63
4 6 1 1 0 UP153 1 11.32
4 6 2 1 0 UP218 1  7.80
4 6 3 1 0 UP290 1  8.44
4 6 4 1 0 UP266 1  9.50
4 6 5 1 0 UP115 1  8.57
4 6 6 1 0 UP301 1 10.81
4 7 1 1 0 UP150 1  6.89
4 7 2 1 0 UP300 1  9.69
4 7 3 1 0 UP282 1 10.53
4 7 4 1 0 UP103 1  7.47
4 7 5 1 0 UP238 1  8.64
4 7 6 1 0 UP151 1  8.73
4 8 1 1 0 UP225 1  9.04
4 8 2 1 0 UP111 1  7.04
4 8 3 1 0 UP283 1  8.62
4 8 4 1 0 UP154 1  6.98
4 8 5 1 0 UP124 1  8.72
4 8 6 1 0 UP140 1  9.44
4 9 1 1 0 UP138 1  8.86
4 9 2 1 0 UP252 1  9.33
4 9 3 1 0 UP110 1  7.52
4 9 4 1 0 UP171 1  9.40
4 9 5 1 0 UP233 1  7.81
4 9 6 1 0 UP156 1  8.80
4 10 1 1 0 UP317 1  8.48
4 10 2 1 0 UP223 1  8.82
4 10 3 1 0 UP316 1  7.44
4 10 4 1 0 UP254 1 10.32
4 10 5 1 0 UP220 1  8.14
4 10 6 1 0 UP287 1  9.51
5 1 1 1 0 UP165 1  9.15
5 1 2 1 0 UP225 1  8.64
5 1 3 1 0 UP220 1  8.67
5 1 4 1 0 UP252 1  8.08
5 1 5 1 0 UP115 1  7.60
5 1 6 1 0 UP300 1 10.31
5 2 1 1 0 UP153 1 11.73
5 2 2 1 0 UP251 1  9.94
5 2 3 1 0 UP151 1  8.76
5 2 4 1 0 UP274 1  9.11
5 2 5 1 0 UP149 1 10.01
5 2 6 1 0 UP276 1 11.03
5 3 1 0 U6 0 0  6.03
5 3 2 1 0 UP223 1  9.70
5 3 3 1 0 UP233 1  9.00
5 3 4 1 0 UP123 1  9.44
5 3 5 1 0 UP301 1  9.70
5 3 6 0 SSOseed 0 0  6.63
5 4 1 1 0 UP283 1 10.01
5 4 2 1 0 UP238 1  8.63
5 4 3 1 0 UP218 1  9.22
5 4 4 1 0 UP297 1  9.04
5 4 5 1 0 UP227 1  8.85
5 4 6 1 0 UP171 1  9.99
5 5 1 1 0 UP291 1 10.06
5 5 2 1 0 UP240 1  7.93
5 5 3 1 0 UP266 1  9.68
5 5 4 1 0 UP190 1  9.03
5 5 5 1 0 UP110 1  6.58
5 5 6 1 0 UP192 1  8.71
5 6 1 1 0 UP282 1  9.73
5 6 2 1 0 UP254 1 10.13
5 6 3 1 0 UP111 1  8.41
5 6 4 1 0 UP188 1  9.26
5 6 5 1 0 UP138 1  8.83
5 6 6 1 0 UP287 1 10.31
5 7 1 1 0 UP219 1  9.01
5 7 2 1 0 UP147 1 10.90
5 7 3 1 0 UP150 1  9.79
5 7 4 1 0 UP154 1  8.24
5 7 5 1 0 UP140 1  8.46
5 7 6 1 0 UP159 1 11.26
5 8 1 1 0 UP144 1  9.53
5 8 2 1 0 UP290 1  8.64
5 8 3 1 0 UP187 1  7.99
5 8 4 1 0 UP278 1  8.04
5 8 5 1 0 UP316 1  8.66
5 8 6 1 0 UP262 1  8.05
5 9 1 1 0 UP124 1  8.81
5 9 2 0 UG323 0 0  8.34
5 9 3 1 0 UP296 1  8.53
5 9 4 1 0 UP169 1  8.07
5 9 5 1 0 UP103 1  7.49
5 9 6 1 0 UP288 1  9.34
5 10 1 1 0 UP122 1  9.11
5 10 2 1 0 UP317 1  9.24
5 10 3 0 PN14 0 0  8.93
5 10 4 1 0 UP156 1  8.72
5 10 5 1 0 UP135 1  9.12
5 10 6 1 0 UP127 1  9.56
;
run;

%macro fit(ddfm, tag);
ods output Tests1 = t1_&tag Tests3 = t3_&tag CovParms = cp_&tag;
proc mixed data=ex82 method=reml;
  class repl row column contcompf std clone;
  model dbh = repl column contcompf contcompf*std / htype=1,3 ddfm=&ddfm;
  random repl*row repl*column;
  random cv / subject=clone;
run;
proc export data=t1_&tag outfile="type1_ex82_tests1_&tag..csv" dbms=csv replace; run;
proc export data=t3_&tag outfile="type1_ex82_tests3_&tag..csv" dbms=csv replace; run;
proc export data=cp_&tag outfile="type1_ex82_covparms_&tag..csv" dbms=csv replace; run;
%mend;

%fit(kr, kr);
%fit(satterthwaite, satt);

/* Run 3 (added 2026-11-25, parallel to dev/genstat/type1-exam8.2.gen run 2):
   another order of the fixed terms, contcompf first. R (predictmeans, X'X
   hypotheses, KR): contcompf F 68.68 (257.5), repl 2.05 (26.5), column 2.00
   (19.3), contcompf*std 16.86 (207). Genstat's X'V^-1X sequential Wald for
   column would be 10.32 (F 2.06). */
ods output Tests1 = t1_order2;
proc mixed data=ex82 method=reml;
  class repl row column contcompf std clone;
  model dbh = contcompf repl column contcompf*std / htype=1 ddfm=kr;
  random repl*row repl*column;
  random cv / subject=clone;
run;
proc export data=t1_order2 outfile="type1_ex82_tests1_order2.csv" dbms=csv replace; run;

