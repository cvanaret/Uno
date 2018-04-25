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

#   classification OOR2-AN-V-V

param n := 1000;
param nx := 2;
param ny := 4;
param c{i in 1..ny, j in 1..nx} := (i+j)/(2*ny);

var x{1..n-1,1..nx};
var y{1..n,1..ny};

minimize f:
	sum {t in 1..n-1} (sum {j in 1..ny} y[t,j]^2)*((sin(0.5*sum{j in 1..nx} x[t,j]^2))^2 + 1.0) + sum {j in 1..ny} y[n,j]^2;
subject to cons1{t in 1..n-1,j in 1..ny}:
	sin(y[t,j]) + sum {i in 1..nx} c[j,i]*(sin(x[t,i])) -y[t+1,j] = 0;
fix{i in 1..ny} y[1,i] := i/(2*ny);
solve;

