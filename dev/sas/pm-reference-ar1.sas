/* Reference output for predictmeans: Kenward-Roger variants for a covariance
   structure that is not linear in its parameters (AR(1)).
   For AR(1), SAS's default DDFM=KR adds second-derivative terms that depend on
   how the structure is parameterised; KR(FIRSTORDER) leaves them out (as
   predictmeans' "KR-SAS"); KR2 is Kenward and Roger (2009). The results
   pages (PDF) are enough; the CSV files are a bonus.
   Data: dev/sas/ar1.txt (simulated: 2 groups x 6 subjects x 4 times,
   AR(1) rho 0.6, 5 of 48 values deleted at random), inlined below. */

%let out = ~/pm_sas_out;
options dlcreatedir;
libname pmout "&out";
libname pmout clear;

data AR1;
  length Group $ 1;
  input Subj Group $ Time Y @@;
  datalines;
1 A 2 10.29  1 A 3 10.02  1 A 4 10.27  2 A 1 10.34  
2 A 2 9.86  2 A 3 11.68  2 A 4 11.96  3 A 1 10.07  
3 A 3 12.44  3 A 4 11.62  4 A 1 11.58  4 A 2 12.61  
4 A 4 11.74  5 A 1 11.58  5 A 2 11.99  5 A 3 12.17  
5 A 4 12.34  6 A 1 10.22  6 A 2 9.39  6 A 3 9.6  
7 B 1 12.37  7 B 2 12.68  7 B 3 12.56  7 B 4 12.95  
8 B 1 12.1  8 B 2 12.39  8 B 3 13.15  8 B 4 12.98  
9 B 1 10.85  9 B 2 10.36  9 B 3 11.22  9 B 4 14  
10 B 1 13.06  10 B 2 13.09  10 B 3 14.72  10 B 4 13.51  
11 B 1 10.32  11 B 2 9.2  11 B 3 11.97  12 B 1 10.35  
12 B 2 11.33  12 B 3 11.86  12 B 4 12.66  
;
run;

%macro exp(ds, name);
  proc export data=&ds outfile="&out/&name..csv" dbms=csv replace; run;
%mend;

%macro ar1(ddfm=, tag=);
  title "AR1 data, DDFM=&ddfm";
  ods output CovParms=cp Tests3=t3 LSMeans=lsm Diffs=dif;
  proc mixed data=AR1 method=reml;
    class Subj Group Time;
    model Y = Group Time Group*Time / ddfm=&ddfm;
    repeated Time / subject=Subj type=ar(1);
    lsmeans Group Time Group*Time / diff cl;
  run;
  %exp(cp, ar1_&tag._covparms) %exp(t3, ar1_&tag._tests3)
  %exp(lsm, ar1_&tag._lsmeans) %exp(dif, ar1_&tag._diffs)
%mend;

%ar1(ddfm=kr,              tag=KR)
%ar1(ddfm=kr(firstorder),  tag=KRFIRST)
%ar1(ddfm=kr2,             tag=KR2)
%ar1(ddfm=satterthwaite,   tag=SATTERTHWAITE)
%ar1(ddfm=betwithin,       tag=BETWITHIN)
