/* Second reference program for predictmeans: multiplicity-adjusted LS-mean
   differences (Tukey-Kramer) with ADJDFE=ROW and ADJDFE=SOURCE (SAS's
   default), and the confidence limits SAS gives for the differences.
   Data: Yates oats without the Block I x Golden Rain whole plot (68
   observations), where the pairs' Kenward-Roger df differ. Run in SAS
   OnDemand as pm-reference.sas; send the results page (PDF) or the CSV
   files written to &out. */

%let out = ~/pm_sas_out;
options dlcreatedir;
libname pmout "&out";
libname pmout clear;

data Oats;
  length Block $ 3 Variety $ 10 Nitro $ 3;
  input Block $ Variety $ Nitro $ Yield @@;
  datalines;
I   Victory     0.0 111  I   Victory     0.2 130  I   Victory     0.4 157  I   Victory     0.6 174
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

%macro adj(adjdfe=);
  ods output Diffs=dif;
  proc mixed data=Oats method=reml;
    class Block Variety Nitro;
    model Yield = Variety Nitro Variety*Nitro / ddfm=kr;
    random int Variety / subject=Block;
    lsmeans Variety Variety*Nitro / diff cl adjust=tukey adjdfe=&adjdfe;
  run;
  proc export data=dif outfile="&out/oatsmisswp_KR_tukey_&adjdfe..csv" dbms=csv replace; run;
%mend;
%adj(adjdfe=row)
%adj(adjdfe=source)
