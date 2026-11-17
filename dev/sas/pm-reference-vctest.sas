/* Reference likelihood ratio tests of variance components for
   varcomp(test = TRUE) (R/vc-test.R, tests/testthat/test-vc-test.R).
   Each model is fitted in full and reduced form (the -2 Res Log Likelihood of
   each fit gives the LRT), and PROC GLIMMIX's COVTEST gives the LRT with its
   p-value (SAS uses the chi-bar-square mixture on the boundary and says so in
   a note). Run in SAS (or SAS OnDemand / SAS Studio) and send back the results
   page (saved as PDF or HTML). Nothing needs editing.
   Cases: (1) Orthodont, random intercept and slope of age (unstructured):
   drop the term, drop the slope, covariance 0; (2) Oats split plot: Block and
   Block*Variety variances; (3) cbpp (lme4), binomial GLMM with a random
   intercept per herd, Laplace; (4) Orthodont, random intercept + AR(1): the
   correlation 0. */

data Oats;
  /* length first: list input reads character values as 8 characters by
     default, which made 'GoldenRain' 'GoldenRa' (the first run, 2026-11-12,
     then deleted nothing for Oats_MissWP) */
  length Block $ 3 Variety $ 10 Nitro $ 3;
  input Block $ Variety $ Nitro $ Yield @@;
  datalines;
I   Victory     0.0 111  I   Victory     0.2 130  I   Victory     0.4 157  I   Victory     0.6 174
I   GoldenRain  0.0 117  I   GoldenRain  0.2 114  I   GoldenRain  0.4 161  I   GoldenRain  0.6 141
I   Marvellous  0.0 105  I   Marvellous  0.2 140  I   Marvellous  0.4 118  I   Marvellous  0.6 156
II  Victory     0.0  61  II  Victory     0.2  91  II  Victory     0.4  97  II  Victory     0.6 100
II  GoldenRain  0.0  70  II  GoldenRain  0.2 108  II  GoldenRain  0.4 126  II  GoldenRain  0.6 149
II  Marvellous  0.0  96  II  Marvellous  0.2 124  II  Marvellous  0.4 121  II  Marvellous  0.6 144
III Victory     0.0  68  III Victory     0.2  64  III Victory     0.4 112  III Victory     0.6  86
III GoldenRain  0.0  60  III GoldenRain  0.2 102  III GoldenRain  0.4  89  III GoldenRain  0.6  96
III Marvellous  0.0  89  III Marvellous  0.2 129  III Marvellous  0.4 132  III Marvellous  0.6 124
IV  Victory     0.0  74  IV  Victory     0.2  89  IV  Victory     0.4  81  IV  Victory     0.6 122
IV  GoldenRain  0.0  64  IV  GoldenRain  0.2 103  IV  GoldenRain  0.4 132  IV  GoldenRain  0.6 133
IV  Marvellous  0.0  70  IV  Marvellous  0.2  89  IV  Marvellous  0.4 104  IV  Marvellous  0.6 117
V   Victory     0.0  62  V   Victory     0.2  90  V   Victory     0.4 100  V   Victory     0.6 116
V   GoldenRain  0.0  80  V   GoldenRain  0.2  82  V   GoldenRain  0.4  94  V   GoldenRain  0.6 126
V   Marvellous  0.0  63  V   Marvellous  0.2  70  V   Marvellous  0.4 109  V   Marvellous  0.6  99
VI  Victory     0.0  53  VI  Victory     0.2  74  VI  Victory     0.4 118  VI  Victory     0.6 113
VI  GoldenRain  0.0  89  VI  GoldenRain  0.2  82  VI  GoldenRain  0.4  86  VI  GoldenRain  0.6 104
VI  Marvellous  0.0  97  VI  Marvellous  0.2  99  VI  Marvellous  0.4 119  VI  Marvellous  0.6 121
;
run;
data Orthodont;
  length Subject $ 3 Sex $ 6;
  input Subject $ Sex $ age distance @@;
  datalines;
M01 Male 8 26.0   M01 Male 10 25.0   M01 Male 12 29.0   M01 Male 14 31.0
M02 Male 8 21.5   M02 Male 10 22.5   M02 Male 12 23.0   M02 Male 14 26.5
M03 Male 8 23.0   M03 Male 10 22.5   M03 Male 12 24.0   M03 Male 14 27.5
M04 Male 8 25.5   M04 Male 10 27.5   M04 Male 12 26.5   M04 Male 14 27.0
M05 Male 8 20.0   M05 Male 10 23.5   M05 Male 12 22.5   M05 Male 14 26.0
M06 Male 8 24.5   M06 Male 10 25.5   M06 Male 12 27.0   M06 Male 14 28.5
M07 Male 8 22.0   M07 Male 10 22.0   M07 Male 12 24.5   M07 Male 14 26.5
M08 Male 8 24.0   M08 Male 10 21.5   M08 Male 12 24.5   M08 Male 14 25.5
M09 Male 8 23.0   M09 Male 10 20.5   M09 Male 12 31.0   M09 Male 14 26.0
M10 Male 8 27.5   M10 Male 10 28.0   M10 Male 12 31.0   M10 Male 14 31.5
M11 Male 8 23.0   M11 Male 10 23.0   M11 Male 12 23.5   M11 Male 14 25.0
M12 Male 8 21.5   M12 Male 10 23.5   M12 Male 12 24.0   M12 Male 14 28.0
M13 Male 8 17.0   M13 Male 10 24.5   M13 Male 12 26.0   M13 Male 14 29.5
M14 Male 8 22.5   M14 Male 10 25.5   M14 Male 12 25.5   M14 Male 14 26.0
M15 Male 8 23.0   M15 Male 10 24.5   M15 Male 12 26.0   M15 Male 14 30.0
M16 Male 8 22.0   M16 Male 10 21.5   M16 Male 12 23.5   M16 Male 14 25.0
F01 Female 8 21.0   F01 Female 10 20.0   F01 Female 12 21.5   F01 Female 14 23.0
F02 Female 8 21.0   F02 Female 10 21.5   F02 Female 12 24.0   F02 Female 14 25.5
F03 Female 8 20.5   F03 Female 10 24.0   F03 Female 12 24.5   F03 Female 14 26.0
F04 Female 8 23.5   F04 Female 10 24.5   F04 Female 12 25.0   F04 Female 14 26.5
F05 Female 8 21.5   F05 Female 10 23.0   F05 Female 12 22.5   F05 Female 14 23.5
F06 Female 8 20.0   F06 Female 10 21.0   F06 Female 12 21.0   F06 Female 14 22.5
F07 Female 8 21.5   F07 Female 10 22.5   F07 Female 12 23.0   F07 Female 14 25.0
F08 Female 8 23.0   F08 Female 10 23.0   F08 Female 12 23.5   F08 Female 14 24.0
F09 Female 8 20.0   F09 Female 10 21.0   F09 Female 12 22.0   F09 Female 14 21.5
F10 Female 8 16.5   F10 Female 10 19.0   F10 Female 12 19.0   F10 Female 14 19.5
F11 Female 8 24.5   F11 Female 10 25.0   F11 Female 12 28.0   F11 Female 14 28.0
;
run;
data cbpp;
  input herd period incidence size @@;
  datalines;
1 1 2 14   1 2 3 12   1 3 4 9   1 4 0 5
2 1 3 22   2 2 1 18   2 3 1 21   3 1 8 22
3 2 2 16   3 3 0 16   3 4 2 20   4 1 2 10
4 2 0 10   4 3 2 9   4 4 0 6   5 1 5 18
5 2 0 25   5 3 0 24   5 4 1 4   6 1 3 17
6 2 0 17   6 3 0 18   6 4 1 20   7 1 8 16
7 2 1 10   7 3 3 9   7 4 0 5   8 1 12 34
9 1 2 9   9 2 0 6   9 3 0 8   9 4 0 6
10 1 1 22   10 2 1 22   10 3 0 18   10 4 2 22
11 1 0 25   11 2 5 27   11 3 3 22   11 4 1 22
12 1 2 10   12 2 1 8   12 3 0 6   12 4 0 5
13 1 1 21   13 2 2 24   13 3 0 19   13 4 0 23
14 1 11 19   14 2 0 2   14 3 0 3   14 4 0 2
15 1 1 19   15 2 1 15   15 3 1 15   15 4 0 15
;
run;

/* ---- 1. Orthodont: distance = age, random intercept and slope by Subject ---- */
title "vctest 1a: Orthodont, (age | Subject) unstructured (full)";
proc mixed data=Orthodont method=reml;
  class Subject;  model distance = age;
  random int age / subject=Subject type=un;
run;
title "vctest 1b: Orthodont, (1 | Subject): slope and covariance dropped";
proc mixed data=Orthodont method=reml;
  class Subject;  model distance = age;
  random int / subject=Subject;
run;
title "vctest 1c: Orthodont, (age || Subject): covariance 0";
proc mixed data=Orthodont method=reml;
  class Subject;  model distance = age;
  random int age / subject=Subject type=vc;
run;
title "vctest 1d: Orthodont, no random term";
proc mixed data=Orthodont method=reml;
  class Subject;  model distance = age;
run;
title "vctest 1e: Orthodont in GLIMMIX with COVTEST (covparms UN(1,1) UN(2,1) UN(2,2) Residual)";
proc glimmix data=Orthodont;
  class Subject;  model distance = age;
  random int age / subject=Subject type=un;
  covtest 'no random term (zero G)'      zerog;
  covtest 'no slope: UN(2,1) = UN(2,2) = 0' . 0 0 .;
  covtest 'covariance 0: UN(2,1) = 0'     indep;
run;

/* ---- 2. Oats split plot: random Block and Block*Variety ---- */
title "vctest 2a: Oats, random Block + Block*Variety (full)";
proc mixed data=Oats method=reml;
  class Block Variety Nitro;  model Yield = Variety Nitro Variety*Nitro;
  random int Variety / subject=Block;
run;
title "vctest 2b: Oats, Block*Variety dropped";
proc mixed data=Oats method=reml;
  class Block Variety Nitro;  model Yield = Variety Nitro Variety*Nitro;
  random int / subject=Block;
run;
title "vctest 2c: Oats, Block dropped";
proc mixed data=Oats method=reml;
  class Block Variety Nitro;  model Yield = Variety Nitro Variety*Nitro;
  random Block*Variety;
run;
title "vctest 2d: Oats in GLIMMIX with COVTEST (covparms Block, Block*Variety, Residual)";
proc glimmix data=Oats;
  class Block Variety Nitro;  model Yield = Variety Nitro Variety*Nitro;
  random Block Block*Variety;
  covtest 'Block variance 0'         0 . .;
  covtest 'Block*Variety variance 0' . 0 .;
run;

/* ---- 3. cbpp: binomial GLMM, random intercept per herd (Laplace) ---- */
title "vctest 3: cbpp, incidence/size = period, random herd (Laplace), COVTEST herd variance 0";
proc glimmix data=cbpp method=laplace;
  class herd period;
  model incidence/size = period / dist=binomial link=logit;
  random int / subject=herd;
  covtest 'herd variance 0' 0;
run;
title "vctest 3b: cbpp, no random term (logistic regression)";
proc glimmix data=cbpp method=laplace;
  class period;
  model incidence/size = period / dist=binomial link=logit;
run;

/* ---- 4. Orthodont: random intercept + AR(1) within child: correlation 0 ---- */
title "vctest 4a: Orthodont, random intercept + AR(1) residuals (full)";
proc mixed data=Orthodont method=reml;
  class Subject;  model distance = age;
  random int / subject=Subject;
  repeated / subject=Subject type=ar(1);
run;
title "vctest 4b: Orthodont, random intercept, independent residuals (AR(1) correlation 0)";
proc mixed data=Orthodont method=reml;
  class Subject;  model distance = age;
  random int / subject=Subject;
run;
