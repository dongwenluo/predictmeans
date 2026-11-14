/* Reference SEs and confidence limits of variance components for varcomp()
   (tests/testthat/test-varcomp.R), with edge cases: missing plots, a
   random slope with an unstructured covariance, AR(1) with a random intercept,
   an unstructured repeated-measures covariance, AR(1) with missing values, a
   component estimated at zero (Dyestuff2), crossed random effects (Penicillin).
   Run in SAS (or SAS OnDemand / SAS Studio) and send back the results page
   (saved as PDF or HTML). Nothing needs editing. COVTEST prints the SEs (from
   the observed information), CL Satterthwaite limits for variances and Wald
   limits for covariances, ASYCOV the covariance of the estimates. */

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
data Oats_MissWP;  set Oats;  if Block='I' and Variety='GoldenRain' then delete; run;
data Oats_Miss;    set Oats_MissWP; if Block='II' and Variety='Victory' and Nitro='0.6' then delete; run;

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
data Dyestuff;  input Batch $ Yield @@;  datalines;
A 1545   A 1440   A 1440   A 1520   A 1580   B 1540
B 1555   B 1490   B 1560   B 1495   C 1595   C 1550
C 1605   C 1510   C 1560   D 1445   D 1440   D 1595
D 1465   D 1545   E 1595   E 1630   E 1515   E 1635
E 1625   F 1520   F 1455   F 1450   F 1480   F 1445
;
run;
data Dyestuff2; input Batch $ Yield @@;  datalines;
A 7.298   A 3.846   A 2.434   A 9.566   A 7.990   B 5.220
B 6.556   B 0.608   B 11.788   B -0.892   C 0.110   C 10.386
C 13.434   C 5.510   C 8.166   D 2.212   D 4.852   D 7.092
D 9.288   D 4.980   E 0.282   E 9.014   E 4.458   E 9.446
E 7.198   F 1.722   F 4.782   F 8.106   F 0.758   F 3.758
;
run;
data Penicillin; input plate $ sample $ diameter @@;  datalines;
a A 27   a B 23   a C 26   a D 23   a E 23   a F 21
b A 27   b B 23   b C 26   b D 23   b E 23   b F 21
c A 25   c B 21   c C 25   c D 24   c E 24   c F 20
d A 26   d B 23   d C 25   d D 23   d E 23   d F 20
e A 25   e B 22   e C 26   e D 22   e E 23   e F 20
f A 24   f B 22   f C 25   f D 23   f E 22   f F 19
g A 24   g B 20   g C 23   g D 21   g E 22   g F 19
h A 26   h B 22   h C 26   h D 24   h E 24   h F 21
i A 24   i B 21   i C 24   i D 22   i E 22   i F 20
j A 24   j B 21   j C 24   j D 23   j E 22   j F 19
k A 26   k B 23   k C 26   k D 24   k E 24   k F 21
l A 25   l B 22   l C 26   l D 24   l E 24   l F 20
m A 26   m B 24   m C 26   m D 24   m E 25   m F 22
n A 26   n B 23   n C 26   n D 23   n E 23   n F 20
o A 26   o B 23   o C 25   o D 24   o E 24   o F 22
p A 25   p B 22   p C 25   p D 23   p E 23   p F 20
q A 25   q B 21   q C 24   q D 23   q E 23   q F 20
r A 25   r B 22   r C 24   r D 23   r E 23   r F 19
s A 24   s B 21   s C 23   s D 21   s E 21   s F 19
t A 26   t B 23   t C 26   t D 24   t E 24   t F 21
u A 25   u B 21   u C 24   u D 22   u E 22   u F 18
v A 25   v B 22   v C 25   v D 22   v E 22   v F 20
w A 24   w B 21   w C 24   w D 22   w E 24   w F 19
x A 24   x B 21   x C 24   x D 22   x E 21   x F 18
;
run;
data AR1; input Subj Group $ Time Y @@;  datalines;
1 A 2 10.29   1 A 3 10.02   1 A 4 10.27   2 A 1 10.34   2 A 2 9.86   2 A 3 11.68
2 A 4 11.96   3 A 1 10.07   3 A 3 12.44   3 A 4 11.62   4 A 1 11.58   4 A 2 12.61
4 A 4 11.74   5 A 1 11.58   5 A 2 11.99   5 A 3 12.17   5 A 4 12.34   6 A 1 10.22
6 A 2 9.39   6 A 3 9.60   7 B 1 12.37   7 B 2 12.68   7 B 3 12.56   7 B 4 12.95
8 B 1 12.10   8 B 2 12.39   8 B 3 13.15   8 B 4 12.98   9 B 1 10.85   9 B 2 10.36
9 B 3 11.22   9 B 4 14.00   10 B 1 13.06   10 B 2 13.09   10 B 3 14.72   10 B 4 13.51
11 B 1 10.32   11 B 2 9.20   11 B 3 11.97   12 B 1 10.35   12 B 2 11.33   12 B 3 11.86
12 B 4 12.66
;
run;

%macro oats(data);
  title "varcomp: &data (random Block, Block*Variety)";
  proc mixed data=&data method=reml covtest cl asycov;
    class Block Variety Nitro;
    model Yield = Variety Nitro Variety*Nitro;
    random int Variety / subject=Block;
  run;
%mend;
%oats(Oats) %oats(Oats_MissWP) %oats(Oats_Miss)

title "varcomp: Orthodont, random intercept and slope (unstructured)";
proc mixed data=Orthodont method=reml covtest cl asycov;
  class Subject;
  model distance = age;
  random int age / subject=Subject type=un;
run;

title "varcomp: Orthodont, random intercept and AR(1)";
proc mixed data=Orthodont method=reml covtest cl asycov;
  class Subject;
  model distance = age;
  random int / subject=Subject;
  repeated / subject=Subject type=ar(1);
run;

title "varcomp: Orthodont, unstructured repeated measures";
data Orthodont2; set Orthodont; agef = age; run;
proc mixed data=Orthodont2 method=reml covtest cl asycov;
  class Subject agef;
  model distance = agef;
  repeated agef / subject=Subject type=un;
run;

title "varcomp: AR(1) repeated measures with missing values";
proc mixed data=AR1 method=reml covtest cl asycov;
  class Subj Group Time;
  model Y = Group Time Group*Time;
  repeated Time / subject=Subj type=ar(1);
run;

title "varcomp: Dyestuff (random Batch)";
proc mixed data=Dyestuff method=reml covtest cl asycov;
  class Batch;  model Yield = ;  random Batch;
run;

title "varcomp: Dyestuff2 (Batch variance estimated at zero)";
proc mixed data=Dyestuff2 method=reml covtest cl asycov;
  class Batch;  model Yield = ;  random Batch;
run;

title "varcomp: Penicillin (crossed: plate and sample)";
proc mixed data=Penicillin method=reml covtest cl asycov;
  class plate sample;  model diameter = ;  random plate sample;
run;
