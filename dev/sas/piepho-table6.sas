/* Piepho (2023), Biometrical Journal 65(7) e2200290, Table 6:
   lamb-level R2 (%) of the cumulative probit (proportional odds) model with
   a random sire effect, for all 16 subsets of the covariates yr b1 b2 b3.

   For each subset this prints the parts of the R2:
     theta_Xb  = variance of the fixed effects, bias-corrected
               = [b'X'PXb - tr(X'PX covB)] / (n - 1)
     theta_Zu  = sire variance x (n - sum n_i^2 / n) / (n - 1)
     theta_h   = 1 (latent residual variance of the probit link)
   computed twice:
     "match"   X = the covariates of this subset (the correct pairing);
     "first_k" X = the first k columns of (yr b1 b2 b3), k = number of
               covariates in the subset, as when the design of the full model
               is reused for every subset.
   Hypothesis (from R, ordinal::clmm with nAGQ = 10): the published Table 6
   equals "first_k" in every row; "match" differs whenever the subset is not
   yr, yr b1, yr b1 b2 or yr b1 b2 b3.

   Output: printed table and CSV file piepho-table6-results.csv (set the path
   in the macro variable outcsv). Run in SAS 9.4 (GLIMMIX, IML). */

%let outcsv = piepho-table6-results.csv;

data foot_mv;
   input yr b1 b2 b3 k1 k2 k3;
   sire = _n_;
   datalines;
  1  1  0  0  52 25  0
  1  1  0  0  49 17  1
  1  1  0  0  50 13  1
  1  1  0  0  42  9  0
  1  1  0  0  74 15  0
  1  1  0  0  54  8  0
  1  1  0  0  96 12  0
  1 -1  1  0  57 52  9
  1 -1  1  0  55 27  5
  1 -1  1  0  70 36  4
  1 -1  1  0  70 37  3
  1 -1  1  0  82 21  1
  1 -1  1  0  75 19  0
  1 -1 -1  0  17 12 10
  1 -1 -1  0  13 23  3
  1 -1 -1  0  21 17  3
 -1  0  0  1  37 41 23
 -1  0  0  1  47 24 12
 -1  0  0  1  46 25  9
 -1  0  0  1  79 32 11
 -1  0  0  1  50 23  5
 -1  0  0  1  63 18  8
 -1  0  0 -1  30 20  9
 -1  0  0 -1  31 33  3
 -1  0  0 -1  28 18  4
 -1  0  0 -1  42 27  4
 -1  0  0 -1  35 22  2
 -1  0  0 -1  33 18  3
 -1  0  0 -1  35 17  4
 -1  0  0 -1  26 13  2
 -1  0  0 -1  37 15  2
 -1  0  0 -1  36 14  1
 -1  0  0 -1  63 20  3
 -1  0  0 -1  41  8  1
;

/* one row per sire x category, with the count of lambs */
data footshape;
   set foot_mv;
   array k{3};
   n = k1 + k2 + k3;
   do Shape = 1 to 3;
      count = k{Shape};
      output;
   end;
   drop k1-k3;
run;

/* one row per sire with its number of lambs: the lamb-level design is these
   rows weighted by n (all lambs of a sire share the covariates) */
data sires;
   set foot_mv;
   n = k1 + k2 + k3;
   keep sire yr b1 b2 b3 n;
run;

data results; length subset $20; stop; run;

%macro r2sub(vars, k);
   %local subset;
   %let subset = %sysfunc(ifc(%length(&vars) = 0, none, &vars));

   ods exclude all;
   proc glimmix data=footshape method=quad;
      ods output ParameterEstimates=pe CovB=cb CovParms=cp FitStatistics=fs;
      class sire;
      model Shape = &vars / s link=cumprobit dist=multinomial covb;
      random int / sub=sire;
      freq count;
   run;
   ods exclude none;

   proc iml;
      use sires; read all var {yr b1 b2 b3} into Xall; read all var {n} into w; close sires;
      use cp; read all var {Estimate} into s2; close cp;
      use fs; read all var {Value} into fit; close fs;
      N = sum(w);
      /* theta_Zu: sire variance x (N - sum n_i^2 / N) / (N - 1) */
      tZu = s2[1] # (N - ssq(w) / N) / (N - 1);
      %if &k = 0 %then %do;
         tXb_m = 0; tXb_f = 0;
      %end;
      %else %do;
         /* the k slopes and their covariance: rows/columns after the two thresholds */
         use pe; read all var {Estimate} into est; close pe;
         use cb; read all var ("Col1":"Col%eval(&k + 2)") into V; close cb;
         b  = est[3:%eval(&k + 2)];
         Vb = V[3:%eval(&k + 2), 3:%eval(&k + 2)];
         /* weighted centring = P X at lamb level */
         start txb(X, b, Vb, w, N);
            m  = (w` * X) / N;
            Xc = X - repeat(m, nrow(X), 1);
            A  = Xc` * (Xc # w);                     /* X'PX at lamb level */
            return( (b` * A * b - trace(A * Vb)) / (N - 1) );
         finish;
         use sires; read all var {&vars} into Xm; close sires;
         Xf = Xall[, 1:&k];
         tXb_m = txb(Xm, b, Vb, w, N);
         tXb_f = txb(Xf, b, Vb, w, N);
      %end;
      tot_m = tXb_m + tZu + 1;
      tot_f = tXb_f + tZu + 1;
      subset = "&subset";
      theta_Zu = tZu;
      theta_Xb_match = tXb_m;   theta_Xb_first_k = tXb_f;
      Fixed_match   = 100 * tXb_m / tot_m;  Random_match   = 100 * tZu / tot_m;
      Total_match   = Fixed_match + Random_match;
      Fixed_first_k = 100 * tXb_f / tot_f;  Random_first_k = 100 * tZu / tot_f;
      Total_first_k = Fixed_first_k + Random_first_k;
      sire_var = s2[1];
      AIC = fit[2];   /* -2 log L, AIC, ... in FitStatistics order */
      create one var {subset sire_var theta_Zu theta_Xb_match theta_Xb_first_k
                      Fixed_match Random_match Total_match
                      Fixed_first_k Random_first_k Total_first_k AIC};
      append; close one;
   quit;

   data results; set results one; run;
%mend;

/* the 16 subsets, in the order of Table 6 */
%r2sub(, 0);
%r2sub(yr, 1);
%r2sub(b1, 1);
%r2sub(b2, 1);
%r2sub(b3, 1);
%r2sub(yr b1, 2);
%r2sub(yr b2, 2);
%r2sub(yr b3, 2);
%r2sub(b1 b2, 2);
%r2sub(b1 b3, 2);
%r2sub(b2 b3, 2);
%r2sub(yr b1 b2, 3);
%r2sub(yr b1 b3, 3);
%r2sub(yr b2 b3, 3);
%r2sub(b1 b2 b3, 3);
%r2sub(yr b1 b2 b3, 4);

proc print data=results noobs;
   format sire_var theta_Zu theta_Xb_match theta_Xb_first_k 10.6
          Fixed_match Random_match Total_match Fixed_first_k Random_first_k Total_first_k 8.2 AIC 10.2;
run;

proc export data=results outfile="&outcsv" dbms=csv replace; run;
