# AMPL Model by Hande Y. Benson
#
# Copyright (C) 2001 Princeton University
# All Rights Reserved
#
# Permission to use, copy, modify, and distribute this software and
# its documentation for any purpose and without fee is hereby
# granted, provided that the above copyright notice appear in all
# copies and that the copyright notice and this
# permission notice appear in all supporting documentation.                     

#   Source: p. 109 of
#   P.A. Steenbrink,
#   "Optimization of Transport Networks",
#   Wiley, 1974.

#   SIF input: Ph. Toint, June 1990.

#   classification QNR2-AY-432-108

param narcs := 18;
param ntrips := 12;

param cost{1..narcs};

param a{i in 1..narcs} := if (i==4) then 0.0 else cost[i]*0.01;
param b{i in 1..narcs} := cost[i]*0.00001;

param rhs{1..ntrips,1..5} default 0.0;

var d{1..ntrips,1..narcs} >= 0.0, := 0.0;
var r{1..ntrips,1..narcs} >= 0.0, := 0.0;

minimize f:
	sum {i in 1..narcs} (a[i]*sum {j in 1..ntrips} d[j,i] + b[i]*(sum {j in 1..ntrips} d[j,i])^2 +
	a[i]*sum {j in 1..ntrips} r[j,i] + b[i]*(sum {j in 1..ntrips} r[j,i])^2);
subject to cons1{i in 1..ntrips}:
	r[i,1]+r[i,2]+r[i,3]-d[i,1]-d[i,2]-d[i,3] = rhs[i,1];
subject to cons2{i in 1..ntrips}:
	-d[i,4]+r[i,4]-d[i,5]+r[i,5]-d[i,6]+r[i,6]+d[i,1]-r[i,1] = rhs[i,2];
subject to cons3{i in 1..ntrips}:
	-d[i,7]+r[i,7]-d[i,8]+r[i,8]-d[i,9]+r[i,9]+d[i,2]-r[i,2] = rhs[i,3];
subject to cons4{i in 1..ntrips}:
	-d[i,10]+r[i,10]-d[i,11]+r[i,11]-d[i,12]+r[i,12]+d[i,4]-r[i,4] = rhs[i,4];
subject to cons5{i in 1..ntrips}:
	-d[i,13]+r[i,13]-d[i,14]+r[i,14]-d[i,15]+r[i,15]+d[i,7]-r[i,7] = rhs[i,5];
subject to cons6{i in 1..ntrips}:
	-d[i,16]+r[i,16]+d[i,10]-r[i,10]+d[i,13]-r[i,13] = 0.0;
subject to cons7{i in 1..ntrips}:
	-d[i,17]+r[i,17]+d[i,3]-r[i,3]+d[i,5]-r[i,5]+d[i,8]-r[i,8] = 0.0;
subject to cons8{i in 1..ntrips}:
	-d[i,18]+r[i,18]+d[i,6]-r[i,6]+d[i,9]-r[i,9]+d[i,11]-r[i,11]+d[i,14]-r[i,14]+d[i,17]-r[i,17] = 0.0;
subject to cons9{i in 1..ntrips}:
	d[i,12]-r[i,12]+d[i,15]-r[i,15]+d[i,16]-r[i,16]+d[i,18]-r[i,18] = 0.0;

data;

param cost :=
1               35.0
2               40.0
3               30.0
4               100.0
5               15.0
6               55.0
7               100.0
8               25.0
9               60.0
10              35.0
11              55.0
12              15.0
13              40.0
14              60.0
15              25.0
16              30.0
17              50.0
18              50.0
;

param rhs :=
1 2 -2000
1 3 2000
2 2 -2000
2 4 2000
3 2 -1000
3 5 1000
4 3 -1000
4 4 1000
5 3 -2000
5 5 2000
6 4 -1000
6 5 1000
7 3 -200
7 2 200
8 4 -200
8 2 200
9 5 -100
9 2 100
10 4 -100
10 3 100
11 5 -200
11 3 200
12 5 -100
12 4 100
;

solve;
