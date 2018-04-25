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

#   Source: a modification of an example in
#   L.G. Willumsen,
#   "Origin-Destination Matrix: static estimation"
#   in "Concise Encyclopedia of Traffic and Transportation Systems"
#   (M. Papageorgiou, ed.), Pergamon Press, 1991.

#   M. Bierlaire, private communication, 1991.

#   SIF input: Ph Toint, Dec 1991.

#   classification OLR2-MN-10-6

set KNOWN;
set UNKNOWN;
param ARCS;
param TC{1..ARCS};
param QLT{1..ARCS};
param P{KNOWN union UNKNOWN, 1..ARCS};
param APV{KNOWN} default 0;
param MU{KNOWN} default 0;
param GAMMA;
param ENTROP;

var T{j in KNOWN union UNKNOWN} >= 0.1, := if (j in KNOWN) then
				APV[j] else 1.0;
var F{i in 1..ARCS} >= 0.1, := TC[i];

minimize f:
	sum {j in KNOWN} MU[j]*((T[j]*log(T[j]/APV[j]))-T[j])
	+ sum {j in UNKNOWN} ENTROP*((T[j]*log(T[j]))-T[j])
	+ sum {i in 1..ARCS} (QLT[i]/GAMMA)*((F[i]*log(F[i]/TC[i]))
	-F[i])
;
subject to cons1{i in 1..ARCS}:
	sum {j in KNOWN union UNKNOWN} P[j,i]*T[j] - F[i] = 0;

data;
set KNOWN := 13 14 23;
set UNKNOWN = 24;
param ARCS := 6;
param:
	TC	QLT:=
1	100	1
2	500	1
3	400	1
4	1100	1
5	600	1
6	700	1;

param P:
	1	2	3	4	5	6:=
13	1	0	0	0	0	0
14	0	1	0	1	0	0
23	0	0	1	1	1	0
24	0	0	0	1	1	1;

param:
	APV	MU:=
13	90	0.5
14	450	0.5
23	360	0.5;

param GAMMA := 1.5;
param ENTROP := 0.2;

solve;
display f;
display T, F;
