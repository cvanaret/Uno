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

#   Source:
#   L.C.W. Dixon, personnal communication, Jan 1991.

#   SIF input: Ph. Toint, Feb 1991.

#   classification SOR2-AN-10-5

param x_init{1..10};
var x{i in 1..10} := x_init[i];

minimize f:
	sum {i in 1..7} 100*(x[i+1]-x[i]^2)^2
	+ sum {i in 1..7} (x[i]-1)^2
	+ sum {i in 1..7} 90*(x[i+3]-x[i+2]^2)^2
	+ sum {i in 1..7} (x[i+2]-1)^2
	+ sum {i in 1..7} 10.1*(x[i+1]-1)^2
	+ sum {i in 1..7} 10.1*(x[i+3]-1)^2
	+ sum {i in 1..7} 19.8*(x[i+1]-1)*(x[i+3]-1)
;

subject to cons1{i in 2..10 by 2}:
	(prod {j in 1..i} x[j]) - 1 = 0;

data;
param x_init:=
1	-2
2	-0.5
3	3
4	0.33333
5	-4
6	-0.25
7	5
8	0.2
9	-6
10	-0.16667;

solve;
display f;
display x;
