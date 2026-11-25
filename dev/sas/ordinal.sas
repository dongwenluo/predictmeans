/* predictmeans: SAS reference for ordinal models (B1-B3a), 2026-11-25.
   Run the whole file; save the HTML results as one PDF and upload it.
   Every predictmeans number for these models (latent means, pairs, class
   probabilities, mean class, their SEs) is a function of the parameter
   estimates and their covariance, so the tables are printed to 16 digits
   (format best16.) and compared in R (dev/sas/compare-ordinal.R).

   Data: ordinal::wine (Randall 1989; 72 bitterness ratings, 1-5, of wine
   from 2 temperatures x 2 contact levels, 9 judges, 8 bottles), and sim6:
   60 ratings (1-3) in 6 clusters, treatment between clusters (the
   predictmeans test of the clmm df).

   Checks:
   1. Cumulative logit, Newton (observed information): = ordinal::clm.
   2. The same with Fisher scoring (expected information): = VGAM::vglm
      (predictmeans SEs of vglm differ from clm by up to 2%; this shows why).
   3. Probit and cloglog cumulative links.
   4. Adjacent-category logit (LINK=ALOGIT): = VGAM::vglm(acat(parallel = TRUE)).
   5. GLIMMIX cumulative logit with a random judge intercept, Laplace: =
      ordinal::clmm; its containment df for temp and contact (within judges).
   6. GLIMMIX on sim6: containment df for a treatment between 6 clusters
      (predictmeans' inner-outer df: 4).
   7. Predicted class probabilities of the four wine cells (model 1). */

options nodate nonumber;
ods html;

data wine;
  input rating temp $ contact $ bottle judge response;
  datalines;
2 cold no 1 1 36
3 cold no 2 1 48
3 cold yes 3 1 47
4 cold yes 4 1 67
4 warm no 5 1 77
4 warm no 6 1 60
5 warm yes 7 1 83
5 warm yes 8 1 90
1 cold no 1 2 17
2 cold no 2 2 22
1 cold yes 3 2 14
3 cold yes 4 2 50
2 warm no 5 2 30
3 warm no 6 2 51
5 warm yes 7 2 90
4 warm yes 8 2 70
2 cold no 1 3 36
3 cold no 2 3 50
3 cold yes 3 3 42
2 cold yes 4 3 23
5 warm no 5 3 80
5 warm no 6 3 81
4 warm yes 7 3 73
4 warm yes 8 3 62
3 cold no 1 4 46
2 cold no 2 4 27
3 cold yes 3 4 48
2 cold yes 4 4 32
3 warm no 5 4 57
2 warm no 6 4 37
5 warm yes 7 4 84
3 warm yes 8 4 58
2 cold no 1 5 26
3 cold no 2 5 45
4 cold yes 3 5 61
3 cold yes 4 5 41
3 warm no 5 5 48
3 warm no 6 5 41
3 warm yes 7 5 58
3 warm yes 8 5 55
3 cold no 1 6 46
2 cold no 2 6 30
3 cold yes 3 6 54
2 cold yes 4 6 37
2 warm no 5 6 32
4 warm no 6 6 60
5 warm yes 7 6 88
4 warm yes 8 6 73
1 cold no 1 7 13
1 cold no 2 7 19
2 cold yes 3 7 31
2 cold yes 4 7 29
2 warm no 5 7 22
3 warm no 6 7 43
2 warm yes 7 7 32
3 warm yes 8 7 49
2 cold no 1 8 25
2 cold no 2 8 32
2 cold yes 3 8 39
3 cold yes 4 8 40
3 warm no 5 8 51
3 warm no 6 8 45
3 warm yes 7 8 42
4 warm yes 8 8 67
1 cold no 1 9 12
2 cold no 2 9 29
3 cold yes 3 9 47
2 cold yes 4 9 28
3 warm no 5 9 47
2 warm no 6 9 38
4 warm yes 7 9 72
4 warm yes 8 9 65
;
run;

data sim6;
  input g trt $ y;
  datalines;
1 A 2
1 A 2
1 A 3
1 A 3
1 A 1
1 A 2
1 A 3
1 A 1
1 A 1
1 A 1
2 A 1
2 A 1
2 A 1
2 A 3
2 A 2
2 A 3
2 A 2
2 A 3
2 A 2
2 A 2
3 A 1
3 A 2
3 A 2
3 A 2
3 A 3
3 A 1
3 A 3
3 A 1
3 A 2
3 A 3
4 B 1
4 B 2
4 B 1
4 B 1
4 B 1
4 B 1
4 B 1
4 B 2
4 B 1
4 B 1
5 B 3
5 B 3
5 B 3
5 B 3
5 B 2
5 B 1
5 B 1
5 B 2
5 B 3
5 B 1
6 B 3
6 B 1
6 B 3
6 B 2
6 B 3
6 B 2
6 B 2
6 B 1
6 B 3
6 B 3
;
run;

%macro show(title);
  title "&title";
  proc print data=pe noobs; format _numeric_ best16.; run;
  proc print data=cb noobs; format _numeric_ best16.; run;
%mend;

/* 1. cumulative logit, Newton-Raphson (observed information): R clm() */
ods output ParameterEstimates=pe CovB=cb;
proc logistic data=wine;
  class temp(ref='cold') contact(ref='no') / param=ref;
  model rating = temp contact temp*contact / link=clogit technique=newton covb;
  output out=pp predprobs=individual;
run;
%show(1. Cumulative logit, Newton (observed information));

/* 7. predicted class probabilities of the four cells (model 1) */
proc sort data=pp nodupkey out=cells; by temp contact; run;
title "7. Class probabilities of the four cells (model 1)";
proc print data=cells noobs; var temp contact IP_1 IP_2 IP_3 IP_4 IP_5; format IP_: best16.; run;

/* 2. the same model, Fisher scoring (expected information): R vglm() */
ods output ParameterEstimates=pe CovB=cb;
proc logistic data=wine;
  class temp(ref='cold') contact(ref='no') / param=ref;
  model rating = temp contact temp*contact / link=clogit technique=fisher covb;
run;
%show(2. Cumulative logit, Fisher scoring (expected information));

/* 3. probit and complementary log-log cumulative links (Newton) */
ods output ParameterEstimates=pe CovB=cb;
proc logistic data=wine;
  class temp(ref='cold') contact(ref='no') / param=ref;
  model rating = temp contact temp*contact / link=probit technique=newton covb;
run;
%show(3a. Cumulative probit, Newton);

ods output ParameterEstimates=pe CovB=cb;
proc logistic data=wine;
  class temp(ref='cold') contact(ref='no') / param=ref;
  model rating = temp contact temp*contact / link=cloglog technique=newton covb;
run;
%show(3b. Cumulative complementary log-log, Newton);

/* 4. adjacent-category logit, proportional (Newton): R vglm(acat(parallel = TRUE)) */
ods output ParameterEstimates=pe CovB=cb;
proc logistic data=wine;
  class temp(ref='cold') contact(ref='no') / param=ref;
  model rating = temp contact temp*contact / link=alogit technique=newton covb;
run;
%show(4. Adjacent-category logit, Newton);

/* 5. cumulative logit with a random judge intercept, Laplace: R clmm() */
ods output ParameterEstimates=pe CovB=cb CovParms=cp Tests3=t3;
proc glimmix data=wine method=laplace;
  class temp(ref='cold') contact(ref='no') judge;
  model rating = temp contact temp*contact / dist=multinomial link=cumlogit solution covb;
  random intercept / subject=judge;
run;
%show(5. GLIMMIX cumulative logit, random judge intercept, Laplace);
proc print data=cp noobs; format _numeric_ best16.; run;
proc print data=t3 noobs; run;

/* 6. sim6: treatment between 6 clusters; the containment df */
ods output ParameterEstimates=pe CovB=cb CovParms=cp Tests3=t3;
proc glimmix data=sim6 method=laplace;
  class trt(ref='A') g;
  model y = trt / dist=multinomial link=cumlogit solution covb;
  random intercept / subject=g;
run;
%show(6. GLIMMIX sim6, treatment between 6 clusters, Laplace);
proc print data=cp noobs; format _numeric_ best16.; run;
proc print data=t3 noobs; run;

ods html close;
