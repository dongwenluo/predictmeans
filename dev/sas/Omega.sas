%macro omega(beta=, V=, XZ=, G=G_empty, covB=, units='no', outomega=outomega);

data G_empty;
justadummy=0;
run;

data &beta;
set &beta;
if effect in ('Scale', 'Residual') then delete;
run;

/*get length of beta*/
data _NULL_;
	if 0 then set &beta nobs=n_beta;
	call symputx('nrows_beta',n_beta);
	stop;
run;

/*get length of G*/
data _NULL_;
	if 0 then set &G nobs=n_G;
	call symputx('nrows_G',n_G);
	stop;
run;

/*get length of XZ*/
data _NULL_;
	if 0 then set &XZ nobs=n_X;
	call symputx('nrows_X',n_X);
	stop;
run;

/*clean up all input files*/
data &beta;
set &beta;
keep estimate;
run;

data &V;
set &V;
keep col1-col&nrows_X;
run;

%if "&G" ne "G_empty" %then %do;
data &G;
set &G;
keep col1-col&nrows_G;
run;
%end;

data X;
set &XZ;
keep _X1-_X&nrows_beta;
run;

%if "&G" ne "G_empty" %then %do;
data Z;
set &XZ;
keep _Z1-_Z&nrows_G;
run;
%end;

data &covB;
set &covB;
keep col1-col&nrows_beta;
run;

proc iml;

use &beta;
read all var _num_ into beta;
close &beta;

use &V;
read all var _num_ into V;
close &V;

%if "&G" ne "G_empty" %then %do;

use &G;
read all var _num_ into G;
close &G;

use Z;
read all var _num_ into Z;
close Z;

%end;

use X;
read all var _num_ into X;
close X;

use &covB;
read all var _num_ into covB;
close &covB;

n=nrow(X);
P=i(n)-j(n,n,1/n);

theta_V    =trace(P*V)/(n-1);
theta_Xbeta=t(beta)*t(X)*P*X*beta/(n-1) - trace(t(X)*P*X*covB)/(n-1);
%if "&G" ne "G_empty" %then %do;
  %if &units='yes' %then %do;
     drop_units=j(&nrows_G,&nrows_G,0);
     drop_units[1:&nrows_G-&nrows_X,1:&nrows_G-&nrows_X]=i(&nrows_G-&nrows_X);
	 theta_Zu   =trace(P*Z*G*drop_units*t(Z))/(n-1); 
  %end;
  %if &units='no' %then %do;
    theta_Zu   =trace(P*Z*G*t(Z))/(n-1); 
  %end;
%end;
%else %do;
  theta_Zu=0;
%end;

Omega_beta  =theta_Xbeta/(theta_Xbeta+theta_V);
Omega_u     =theta_Zu   /(theta_Xbeta+theta_V);
Omega_beta_u=Omega_beta+Omega_u;

print Omega_beta Omega_u Omega_beta_u;

create outomega var {Omega_beta Omega_u Omega_beta_u};
append;
close outomega;

Omega_beta_percent=Omega_beta*100;
Omega_u_percent=Omega_u*100;
Omega_beta_u_percent=Omega_beta_u*100;

print Omega_beta_percent Omega_u_percent Omega_beta_u_percent;

run;
quit;

%mend;
