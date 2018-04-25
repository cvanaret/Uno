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

#   Source: problem P4 in
#   W.W. Hager,
#   "Multiplier Methods for Nonlinear Optimal Control",
#   SIAM J. on Numercal Analysis 27(4): 1061-1080, 1990.

#   SIF input: Ph. Toint, April 1991.

#   classification OLR2-AN-V-V

param n := 5000;
param h := 1/n;

param t{i in 0..n} := i*h;
param z{i in 0..n} := exp(-2*t[i]);
param a{i in 0..1} := -0.5*z[i];
param b{i in 0..1} := a[i]*(t[i]+0.5);
param c{i in 0..1} := a[i]*(t[i]^2+t[i]+0.5);

param scda := (a[1] - a[0])/2;
param scdb := (b[1] - b[0])*n;
param scdc := (c[1] - c[0])*n*n*0.5;

param e := exp(1);
param xx0 := (1+3*e)/(2-2*e);

var x{0..n} := 0.0;
var u{1..n} <= 1.0, := 0.0;

minimize f:
	sum {i in 1..n} (scda*z[i-1]*x[i]^2 + scdb*z[i-1]*x[i]*(x[i-1]-x[i]) + scdc*z[i-1]*(x[i-1]-x[i])^2) +
	sum {i in 1..n} (u[i]^2)*h*0.5;
subject to cons1{i in 1..n}:
	(n-1)*x[i] - n*x[i-1] - exp(t[i])*u[i] = 0;

fix x[0] := xx0;
solve; 

