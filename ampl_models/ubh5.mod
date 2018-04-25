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

#   Source: unscaled problem 5 
#   (ODE = 1, CLS = 2, GRD = 1, MET = T, SEED = 0.) in
#   J.T. Betts and W.P. Huffman,
#   "Sparse Nonlinear Programming Test Problems (Release 1.0)",
#   Boeing Computer services, Seattle, July 1993.

#   SIF input: Ph.L. Toint, October 1993.

#   classification LQR2-MN-V-V

param n := 2000;
param t0 := 0.0;
param tf := 1000.0;

param k := (tf-t0)/n;

var y{i in 1..7, t in 0..n} := 0.0;
var u{i in 1..3, t in 0..n} >= -1.0, <= 1.0, := 0.0;

minimize f:
	y[7,n];
subject to cons1{i in 1..3, t in 1..n}:
	y[i,t] - y[i,t-1] - k*y[i+3,t-1]/2 - k*y[i+3,t]/2 = 0;
subject to cons2{i in 1..3, t in 1..n}:
	y[i+3,t] - y[i+3,t-1] - k*u[i,t-1]/2 - k*u[i,t]/2 = 0;
subject to cons3{t in 1..n}:
	sum {i in 1..3} (-k)*(u[i,t]^2+u[i,t-1]^2)/2 + y[7,t] - y[7,t-1]= 0;

fix y[1,0] := 1000.0;
fix y[2,0] := 1000.0;
fix y[3,0] := 1000.0;
fix y[4,0] := -10.0;
fix y[5,0] := 10.0;
fix y[6,0] := -10.0;
fix y[7,0] := 0.0;
fix {i in 1..6} y[i,n] := 0.0;

solve;
