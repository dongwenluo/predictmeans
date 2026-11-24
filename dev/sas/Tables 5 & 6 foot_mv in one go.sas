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
data footshape_binomial; 
set foot_mv;
n=k1+k2+k3;
m=k1;
output;
run;

/*R2 at sire level*/

proc glimmix data=footshape_binomial method=quad;
ods output ParameterEstimates=my_beta covB=my_covB CovParms=cp;
class Sire;
model m/n = yr b1 b2 b3 / s link=probit dist=binomial covb;
random int / sub=Sire;
output out=eta_data pred(blup noilink)=eta variance(ilink) / derivatives;
run;

data eta_data;
set eta_data;
R_h=Variance_Dep/d_MU_wrt_LINP_**2;
inv_R_h=1/R_h;
run;

data scale;
estimate=1;
run;

data cp2;
set cp scale;
run;

ods exclude all;
proc glimmix data=eta_data outdesign(X Z novar)=my_XZ lognote noprofile;
ods output G=my_G V=my_V;
weight inv_R_h;
class Sire;
model k1=yr b1 b2 b3;
random Sire/V G;
parms / pdata=cp2 hold=1,2;
run;
ods exclude none;

%omega(beta=my_beta, XZ=my_XZ, V=my_V, G=my_G, covB=my_covB);

/*R2 at lamb level*/

data footshape_binomial;
set footshape_binomial;
R_h=1;
inv_R_h=1/R_h;
run;

data scale;
estimate=1;
run;

data cp2;
set cp scale;
run;

ods exclude all;
proc glimmix data=footshape_binomial outdesign(X Z novar)=my_XZ lognote noprofile;
ods output G=my_G V=my_V;
weight inv_R_h;
class Sire;
freq n;
model k1=yr b1 b2 b3;
random Sire/V G;
parms / pdata=cp2 hold=1,2;
run;
ods exclude none;

data my_XZ_lamb;
merge footshape_binomial my_XZ;
do i=1 to n;
  output;
end;
run;

%omega(beta=my_beta, XZ=my_XZ_lamb, V=my_V, G=my_G, covB=my_covB);

/*now do the same with the multinomial ordered probit model*/

data footshape_binomial2; 
set foot_mv;
m=k1; n=k1+k2+k3; Shape=1; output;
m=k2; n=k2+k3;    Shape=2; output;
run;

/*R2 at sire level*/

proc glimmix data=footshape_binomial2 method=quad;
ods output ParameterEstimates=my_beta covB=my_covB CovParms=cp;
class Sire Shape;
model m/n = Shape yr b1 b2 b3 / noint s link=probit dist=binomial covb;
random int / sub=Sire;
output out=eta_data pred(blup noilink)=eta variance(ilink) / derivatives;
run;

data eta_data;
set eta_data;
R_h=Variance_Dep/d_MU_wrt_LINP_**2;
inv_R_h=1/R_h;
if Shape=1; /*!!!*/
run;

data scale;
estimate=1;
run;

data cp2;
set cp scale;
run;

ods exclude all;
proc glimmix data=eta_data outdesign(X Z novar)=my_XZ lognote noprofile;
ods output G=my_G V=my_V;
weight inv_R_h;
class Sire;
model k1=yr b1 b2 b3;
random Sire/V G;
parms / pdata=cp2 hold=1,2;
run;
ods exclude none;

data my_beta;
set my_beta;
if _N_=2 then delete; /*drop the second intercept*/
run;

data my_covB;
set my_covB;
if _N_=2 then delete;
drop col2;
rename col3-col6=col2-col5;
run;

%omega(beta=my_beta, XZ=my_XZ, V=my_V, G=my_G, covB=my_covB);

/*now use the proportional odds model*/

/*R2 at lamb level*/

data footshape; set foot_mv;
array k{3};
n=k1+k2+k3;
do Shape = 1 to 3;
   count = k{Shape};
   output;
end;
drop k:;
run;

proc glimmix data=footshape method=quad;
ods output ParameterEstimates=my_beta covB=my_covB CovParms=cp;
class Sire;
model Shape = yr b1 b2 b3 / s link=cumprobit dist=multinomial covb;
random int / sub=Sire;
freq count;
output out=eta_data pred(blup noilink)=eta variance(ilink) / derivatives;
run;

data footshape;
set footshape;
R_h=1;
inv_R_h=1/R_h;
run;

data scale;
estimate=1;
run;

data cp2;
set cp scale;
run;

ods exclude all;
proc glimmix data=footshape outdesign(X Z novar)=my_XZ lognote noprofile;
ods output G=my_G V=my_V;
weight inv_R_h;
class Sire;
freq n;
model Shape=yr b1 b2 b3;
random Sire/V G;
parms / pdata=cp2 hold=1,2;
run;
ods exclude none;

data my_XZ_lamb;
merge footshape my_XZ;
do i=1 to n;
  output;
end;
run;

data my_beta;
set my_beta;
if _N_=2 then delete; /*drop the second intercept*/
run;

data my_covB;
set my_covB;
if _N_=2 then delete;
drop col2;
rename col3-col6=col2-col5;
run;

%omega(beta=my_beta, XZ=my_XZ_lamb, V=my_V, G=my_G, covB=my_covB);

/*sire level*/

data eta_data;
set eta_data;
mu=probnorm(eta);
Variance_Dep=n*mu*(1-mu);
d_MU_wrt_LINP_=n*d_MU_wrt_LINP_;
R_h=Variance_Dep/d_MU_wrt_LINP_**2;
inv_R_h=1/R_h;
if Shape=1;
if _Level_=1;
run;

ods exclude all;
proc glimmix data=eta_data outdesign(X Z novar)=my_XZ lognote noprofile;
ods output G=my_G V=my_V;
weight inv_R_h;
class Sire;
model Shape=yr b1 b2 b3;
random Sire/V G;
parms / pdata=cp2 hold=1,2;
run;
ods exclude none;

%omega(beta=my_beta, XZ=my_XZ, V=my_V, G=my_G, covB=my_covB);
