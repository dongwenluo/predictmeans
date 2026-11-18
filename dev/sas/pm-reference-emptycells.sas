/* Reference F tests for anova() with an empty cell (Type III on the cells
   with data; R/anova.R pm_type3_map_cells(), tests/testthat/test-anova-empty.R).
   Oats with the Victory x 0.6 plots removed in every block (11 cells). The
   model is written on the cell means (factor NV, one level per observed
   cell), so each hypothesis is a CONTRAST on the cells, with nothing left to
   SAS's choice of Type III/IV basis:
     Variety   GR - Vic and Marv - Vic over Nitro 0, 0.2, 0.4 (the Nitro
               levels all varieties have; 2 df)
     Nitro     0, 0.2, 0.4 against 0.6 over Golden Rain and Marvellous (the
               varieties all Nitro levels have; 3 df)
     N x V     the interaction contrasts of the observed cells (5 df)
   (1) PROC GLM, Block fixed: F tests from the residual mean square; the
   factorial model with SS3 and SS4 for comparison. (2) PROC MIXED, random
   Block and Block*Variety, DDFM=KR (KR(FIRSTORDER) too) and SATTERTHWAITE.
   Run in SAS (or SAS OnDemand / SAS Studio) and send back the results page
   (saved as PDF or HTML). Nothing needs editing.
   Cell order of NV: GoldenRain 0.0 0.2 0.4 0.6, Marvellous 0.0 0.2 0.4 0.6,
   Victory 0.0 0.2 0.4. */

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

data OatsE;
  set Oats;
  if Variety = "Victory" and Nitro = "0.6" then delete;
  length NV $ 14;
  NV = catx("_", Variety, Nitro);
run;

%macro cellcontrasts;
  contrast "Variety (common Nitro)"
    NV 1 1 1 0  0 0 0 0  -1 -1 -1,
    NV 0 0 0 0  1 1 1 0  -1 -1 -1;
  contrast "Nitro (common Variety)"
    NV 1 0 0 -1  1 0 0 -1  0 0 0,
    NV 0 1 0 -1  0 1 0 -1  0 0 0,
    NV 0 0 1 -1  0 0 1 -1  0 0 0;
  contrast "Nitro x Variety (observed cells)"
    NV 1 0 0 -1  -1 0 0 1  0 0 0,
    NV 0 1 0 -1  0 -1 0 1  0 0 0,
    NV 0 0 1 -1  0 0 -1 1  0 0 0,
    NV 1 0 -1 0  0 0 0 0  -1 0 1,
    NV 0 1 -1 0  0 0 0 0  0 -1 1;
%mend;

title "(1) PROC GLM, Block fixed: cell-means contrasts";
proc glm data=OatsE;
  class Block NV;
  model Yield = Block NV / ss3;
  %cellcontrasts;
run; quit;

title "(1b) PROC GLM, factorial model: Type III and IV";
proc glm data=OatsE;
  class Block Variety Nitro;
  model Yield = Block Nitro Variety Nitro*Variety / ss3 ss4;
run; quit;

title "(2) PROC MIXED, DDFM=KR";
proc mixed data=OatsE method=reml;
  class Block Variety NV;
  model Yield = NV / ddfm=kr;
  random Block Block*Variety;
  %cellcontrasts;
run;

title "(2b) PROC MIXED, DDFM=KR(FIRSTORDER)";
proc mixed data=OatsE method=reml;
  class Block Variety NV;
  model Yield = NV / ddfm=kr(firstorder);
  random Block Block*Variety;
  %cellcontrasts;
run;

title "(2c) PROC MIXED, DDFM=SATTERTHWAITE";
proc mixed data=OatsE method=reml;
  class Block Variety NV;
  model Yield = NV / ddfm=satterthwaite;
  random Block Block*Variety;
  %cellcontrasts;
run;
