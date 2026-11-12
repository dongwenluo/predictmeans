/* Reference output for predictmeans df methods.
   Run in SAS (or SAS OnDemand / SAS Studio) and send back the folder of CSV
   files written to &out. Nothing here needs editing except the out path.
   Data: Yates oats (dev/sas/oats.txt) and the repeated-measures trial
   (dev/sas/rep.txt), inlined below. */

%let out = ~/pm_sas_out;   /* folder for the CSV files, in your SAS home (SAS OnDemand: Files > Home) */
options dlcreatedir;       /* the LIBNAME below creates the folder if it does not exist */
libname pmout "&out";
libname pmout clear;

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

data RepeatedTrial;
  input Subj Group $ Time1 Time2 Time3;
  datalines;
1  TrtA  20.1  24.5  28.2
2  TrtA  18.4  22.0  25.9
3  TrtA  22.3  26.1  30.4
4  TrtA  19.0  23.4  27.1
5  TrtA  21.5  25.0  29.0
6  TrtA  17.8  21.2  24.8
7  TrtB  19.5  21.0  22.4
8  TrtB  21.0  22.8  23.9
9  TrtB  18.2  19.5  20.8
10 TrtB  20.4  21.9  23.1
11 TrtB  22.1  23.5  24.6
12 TrtB  19.0  20.1  21.5
;
run;
data RepeatedLong; set RepeatedTrial;
  Time=1; Y=Time1; output; Time=2; Y=Time2; output; Time=3; Y=Time3; output;
  keep Subj Group Time Y;
run;

%macro exp(ds, name);
  proc export data=&ds outfile="&out/&name..csv" dbms=csv replace; run;
%mend;

%macro oats(data=, tag=, ddfm=);
  ods output CovParms=cp Tests3=t3 LSMeans=lsm Diffs=dif SolutionF=sol;
  proc mixed data=&data method=reml;
    class Block Variety Nitro;
    model Yield = Variety Nitro Variety*Nitro / ddfm=&ddfm solution;
    random int Variety / subject=Block;
    lsmeans Variety Nitro Variety*Nitro / diff cl;
  run;
  %exp(cp,  &tag._&ddfm._covparms) %exp(t3, &tag._&ddfm._tests3)
  %exp(lsm, &tag._&ddfm._lsmeans)  %exp(dif, &tag._&ddfm._diffs) %exp(sol, &tag._&ddfm._solution)
%mend;

%macro rep(ddfm=);
  ods output CovParms=cp Tests3=t3 LSMeans=lsm Diffs=dif;
  proc mixed data=RepeatedLong method=reml;
    class Subj Group Time;
    model Y = Group Time Group*Time / ddfm=&ddfm;
    repeated Time / subject=Subj type=un;
    lsmeans Group Time Group*Time / diff cl;
  run;
  %exp(cp, rep_&ddfm._covparms) %exp(t3, rep_&ddfm._tests3)
  %exp(lsm, rep_&ddfm._lsmeans) %exp(dif, rep_&ddfm._diffs)
%mend;

%macro all;
  %let ms = CONTAIN BETWITHIN SATTERTHWAITE KR;
  %do i = 1 %to 4;
    %let m = %scan(&ms, &i);
    %oats(data=Oats,        tag=oats,     ddfm=&m)
    %oats(data=Oats_MissWP, tag=misswp,   ddfm=&m)
    %oats(data=Oats_Miss,   tag=miss,     ddfm=&m)
    %rep(ddfm=&m)
  %end;
%mend;
%all
