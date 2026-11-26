/* Inner-outer df of aggregated binomial GLMMs (2026-11-26): what does GLIMMIX
   count, rows or trials? The data of tests/testthat/test-inner-outer.R
   (io_clusters(): 8 clusters, Trt between clusters, W within, 10 trials a row).
   predictmeans: Trt 6 df; W and Trt*W 470 df (trials), 38 with
   options(predictmeans.inner_outer_count = "rows").
   Please run and upload the output (PDF or listing). */
data io;
  input clus Trt $ W unit r n;
  datalines;
1 a 1 1 7 10
1 a 1 2 4 10
1 a 1 3 7 10
1 a 2 1 5 10
1 a 2 2 5 10
1 a 2 3 8 10
2 a 1 1 3 10
2 a 1 2 3 10
2 a 1 3 4 10
2 a 2 1 5 10
2 a 2 2 3 10
2 a 2 3 5 10
3 a 1 1 5 10
3 a 1 2 5 10
3 a 1 3 5 10
3 a 2 1 4 10
3 a 2 2 7 10
3 a 2 3 5 10
4 a 1 1 4 10
4 a 1 2 7 10
4 a 1 3 4 10
4 a 2 1 5 10
4 a 2 2 4 10
4 a 2 3 1 10
5 b 1 1 3 10
5 b 1 2 7 10
5 b 1 3 6 10
5 b 2 1 7 10
5 b 2 2 4 10
5 b 2 3 8 10
6 b 1 1 7 10
6 b 1 2 7 10
6 b 1 3 6 10
6 b 2 1 5 10
6 b 2 2 4 10
6 b 2 3 6 10
7 b 1 1 4 10
7 b 1 2 4 10
7 b 1 3 4 10
7 b 2 1 2 10
7 b 2 2 1 10
7 b 2 3 1 10
8 b 1 1 5 10
8 b 1 2 4 10
8 b 1 3 4 10
8 b 2 1 1 10
8 b 2 2 5 10
8 b 2 3 6 10
;
run;
title "events/trials, DDFM=CONTAIN (default)";
proc glimmix data=io method=quad;
  class clus Trt W;
  model r/n = Trt W Trt*W / dist=binomial link=logit;
  random intercept / subject=clus;
run;
title "events/trials, DDFM=BETWITHIN";
proc glimmix data=io method=quad;
  class clus Trt W;
  model r/n = Trt W Trt*W / dist=binomial link=logit ddfm=betwithin;
  random intercept / subject=clus;
run;
/* the same data, one row per trial */
data io1;
  set io;
  do k = 1 to n; y = (k <= r); output; end;
run;
title "one row per trial, DDFM=BETWITHIN";
proc glimmix data=io1 method=quad;
  class clus Trt W;
  model y(event = "1") = Trt W Trt*W / dist=binary link=logit ddfm=betwithin;
  random intercept / subject=clus;
run;
