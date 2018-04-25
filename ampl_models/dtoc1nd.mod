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

#   classification OQR2-AN-V-V

param n := 50;
param nx := 5;
param ny := 10;

param mu := 1.0;

param b{i in 1..ny, j in 1..nx} := (i-j)/(nx+ny);
param c{i in 1..ny, j in 1..nx} := (i+j)*mu/(nx+ny);

var x{1..n-1, 1..nx} := 0.0;
var y{1..n, 1..ny} := 0.0;

minimize f:
	sum {t in 1..n-1, i in 1..nx} (x[t,i] + 0.5)^4 +
	sum {t in 1..n, i in 1..ny} (y[t,i] + 0.25)^4;
subject to cons1{t in 1..n-1}:
	sum {k in 0..ny*nx-1} c[(k div nx)+1,k-nx*(k div nx)+1]*y[t,(k div nx)+1]*x[t,k-nx*(k div nx)+1]+0.5*y[t,1] + 0.25*y[t,2] - y[t+1,1] + sum {i in 1..nx} b[1,i]*x[t,i] = 0;
subject to cons2{t in 1..n-1, j in 2..ny-1}:
        sum {k in 0..ny*nx-1} c[(k div nx)+1,k-nx*(k div nx)+1]*y[t,(k div nx)+1]*x[t,k-nx*(k div nx)+1]-y[t+1,j] + 0.5*y[t,j] - 0.25*y[t,j-1] + 0.25*y[t,j+1] + sum {i in 1..nx} b[j,i]*x[t,i] = 0;
subject to cons3{t in 1..n-1}:
        sum {k in 0..ny*nx-1} c[(k div nx)+1,k-nx*(k div nx)+1]*y[t,(k div nx)+1]*x[t,k-nx*(k div nx)+1]+0.5*y[t,ny] - 0.25*y[t,ny-1] - y[t+1,ny] + sum {i in 1..nx} b[ny,i]*x[t,i] = 0;

fix {i in 1..ny} y[1,i] := 0.0;
solve;
