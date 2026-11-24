/* se_ranef() reference: SEs of the predicted random effects (BLUPs) in SAS
   PROC MIXED, to check predictmeans::se_ranef().

   Data: Yates' oats split plot (nlme::Oats; 6 blocks, 3 varieties on whole
   plots, 4 nitrogen rates on subplots), made unbalanced by deleting 7 of the
   72 plots (65 left): the whole plot Block I x Marvellous (4 subplots) and
   the subplots III-Victory-0.4, V-GoldenRain-0.0 and VI-Victory-0.6. With
   missing plots the SEs differ between levels, a sharper test than the
   balanced data (where every level has the same SE).

   Model: yield ~ nitro*Variety + (1 | Block/Variety), REML.

   What to compare (StdErrPred in SolutionR):
     se_ranef(fm)                       type = "prediction": SAS StdErrPred
       (prediction error SE, sqrt Var(u_hat - u), from the inverse of the
        mixed model equations: includes the uncertainty of the fixed effects)
     se_ranef(fm, type = "conditional") lme4 condVar (1.x se_ranef): not in
       SAS; listed for reference

   R values (lme4 REML; variance components Block 284.9621,
   Block:Variety 110.9739, residual 171.5319):
     group         level          estimate  se_pred  se_cond
     Variety:Block GoldenRain:I     -0.3260   8.4076   7.9022
     Variety:Block GoldenRain:II     4.8225   7.9730   7.3173
     Variety:Block GoldenRain:III   -6.2935   8.0025   7.3552
     Variety:Block GoldenRain:IV     6.5881   7.9730   7.3173
     Variety:Block GoldenRain:V      0.1303   8.1482   7.5585
     Variety:Block GoldenRain:VI    -4.9214   8.0025   7.3552
     Variety:Block Marvellous:II     4.9203   8.0688   7.3173
     Variety:Block Marvellous:III   10.9347   8.0946   7.3552
     Variety:Block Marvellous:IV    -8.4610   8.0688   7.3173
     Variety:Block Marvellous:V    -10.3241   8.0938   7.3552
     Variety:Block Marvellous:VI     2.9302   8.0946   7.3552
     Variety:Block Victory:I        11.9333   8.4103   7.9022
     Variety:Block Victory:II       -8.7040   7.9793   7.3173
     Variety:Block Victory:III      -7.9207   8.1594   7.5585
     Variety:Block Victory:IV       -0.0862   7.9793   7.3173
     Variety:Block Victory:V         5.4437   8.0080   7.3552
     Variety:Block Victory:VI       -0.6661   8.1594   7.5585
     Block         I                29.8054   9.9082   7.7830
     Block         II                2.6673   9.1854   6.5927
     Block         III              -8.4212   9.2349   6.6733
     Block         IV               -5.0306   9.1854   6.5927
     Block         V               -12.1974   9.2339   6.6733
     Block         VI               -6.8235   9.2349   6.6733

   R code for the same fit:
     O <- data.frame(nlme::Oats)  # drop the 7 plots above, nitro as factor
     fm <- lmer(yield ~ nitro*Variety + (1|Block/Variety), data = D)
     se_ranef(fm); se_ranef(fm, type = "conditional")

   Run, then send se_ranef_solutionr.csv and se_ranef_covparms.csv. */

data oats;
  length Block $ 3 Variety $ 10 Nitro $ 3;
  input Block $ Variety $ Nitro $ Yield;
  datalines;
I    Victory     0.0  111
I    Victory     0.2  130
I    Victory     0.4  157
I    Victory     0.6  174
I    GoldenRain  0.0  117
I    GoldenRain  0.2  114
I    GoldenRain  0.4  161
I    GoldenRain  0.6  141
II   Victory     0.0   61
II   Victory     0.2   91
II   Victory     0.4   97
II   Victory     0.6  100
II   GoldenRain  0.0   70
II   GoldenRain  0.2  108
II   GoldenRain  0.4  126
II   GoldenRain  0.6  149
II   Marvellous  0.0   96
II   Marvellous  0.2  124
II   Marvellous  0.4  121
II   Marvellous  0.6  144
III  Victory     0.0   68
III  Victory     0.2   64
III  Victory     0.6   86
III  GoldenRain  0.0   60
III  GoldenRain  0.2  102
III  GoldenRain  0.4   89
III  GoldenRain  0.6   96
III  Marvellous  0.0   89
III  Marvellous  0.2  129
III  Marvellous  0.4  132
III  Marvellous  0.6  124
IV   Victory     0.0   74
IV   Victory     0.2   89
IV   Victory     0.4   81
IV   Victory     0.6  122
IV   GoldenRain  0.0   64
IV   GoldenRain  0.2  103
IV   GoldenRain  0.4  132
IV   GoldenRain  0.6  133
IV   Marvellous  0.0   70
IV   Marvellous  0.2   89
IV   Marvellous  0.4  104
IV   Marvellous  0.6  117
V    Victory     0.0   62
V    Victory     0.2   90
V    Victory     0.4  100
V    Victory     0.6  116
V    GoldenRain  0.2   82
V    GoldenRain  0.4   94
V    GoldenRain  0.6  126
V    Marvellous  0.0   63
V    Marvellous  0.2   70
V    Marvellous  0.4  109
V    Marvellous  0.6   99
VI   Victory     0.0   53
VI   Victory     0.2   74
VI   Victory     0.4  118
VI   GoldenRain  0.0   89
VI   GoldenRain  0.2   82
VI   GoldenRain  0.4   86
VI   GoldenRain  0.6  104
VI   Marvellous  0.0   97
VI   Marvellous  0.2   99
VI   Marvellous  0.4  119
VI   Marvellous  0.6  121
;
run;

ods output SolutionR = solr CovParms = cp SolutionF = solf;
proc mixed data=oats method=reml;
  class Block Variety Nitro;
  model Yield = Nitro|Variety / solution;
  random int / subject=Block solution;
  random Variety / subject=Block solution;
run;

proc print data=cp; run;
proc print data=solr; run;

proc export data=solr outfile='se_ranef_solutionr.csv' dbms=csv replace; run;
proc export data=cp   outfile='se_ranef_covparms.csv'  dbms=csv replace; run;
