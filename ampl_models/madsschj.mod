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
#   K. Madsen and H. Schjaer-Jacobsen,
#   "Linearly Constrained Minmax Optimization",
#   Mathematical Programming 14, pp. 208-223, 1978.

#   SIF input: Ph. Toint, August 1993.

#   classification LQR2-AN-V-V

param N:=80;
param M:=2*N-2;

var x{1..N} := 10;
var z := 0;

minimize f:
	z;
subject to cons1:
	z-sum {i in 2..N} x[i]+1-x[1]^2 >= 0;
subject to cons2:
	z-x[1]-sum {i in 3..N} x[i]+1-x[1]^2 >= 0;
subject to cons3:
	z-x[1]-sum {i in 3..N} x[i]+1-2*x[2]^2 >= 0;
subject to cons4{k in 4..M-1 by 2}:
	z-sum {i in 1..(k+2)/2-1} x[i] - sum {i in
	(k+2)/2+1..N} x[i] +1-x[(k+2)/2]^2>= 0;
subject to cons5{k in 4..M-1 by 2}:
	z-sum {i in 1..(k+2)/2-1} x[i] - sum {i in
	(k+2)/2+1..N} x[i] +1-2*x[(k+2)/2]^2>= 0;
subject to cons6:
	z-sum {i in 1..N-1} x[i] +1 - x[N]^2>= 0;

solve;
display f;
display x; display z;
