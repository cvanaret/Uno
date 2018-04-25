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
#   M. Gawande and J. Dunn,
#   "A Projected Newton Method in a Cartesian Product of Balls",
#   JOTA 59(1): 59-69, 1988.

#   SIF input: Ph. Toint, June 1990.

#   classification QQR2-AN-V-V

param n := 10;

param speed := 0.01;
param pen := 0.335;

var x{j in 1..2, i in 0..n+1};
var v{j in 1..2, i in 0..n+1};
var f{j in 1..2, i in 0..n};

minimize obj:
	pen*(v[1,n+1]^2+v[2,n+1]^2) - (x[1,n+1]^2+x[2,n+1]^2);
subject to cons1{i in 1..n+1, j in 1..2}:
	x[j,i] - x[j,i-1] - v[j,i-1]/n - f[j,i-1]/(2*n^2) = 0;
subject to cons2{i in 1..n+1, j in 1..2}:
	v[j,i] - v[j,i-1] - f[j,i-1]/n = 0;
subject to cons3{i in 0..n}:
	f[1,i]^2 + f[2,i]^2 <= 1;

fix x[1,0] := 0.0;
fix x[2,0] := 0.0;
fix v[1,0] := speed;
fix v[2,0] := 0.0;
solve;


