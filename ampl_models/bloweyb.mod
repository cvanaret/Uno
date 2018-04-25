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

#   Source: a simplification of
#   J.F. Blowey and C.M. Elliott,
#   "The Cahn-Hilliard gradient theory for phase separation with 
#   non-smooth free energy Part II: Numerical analysis",
#   European Journal of Applied Mathematics (3) pp 147-179, 1992.

#   SIF input: Nick Gould, August 1996

#   classification QLR2-MN-V-V

param N := 1000;

param A := 0.1;
param B := 0.4;
param C := 0.6;
param D := 0.9;

param INT := N*(1+A+B-C-D);

param v{i in 0..N} := if (0 <= i <= N*A) then 1.0
		else  if (N*A+1 <= i <= N*B) then (1 - (i-N*A)*2/(N*(B-A)))
		else  if (N*B+1 <= i <= N*C) then -1 
		else  if (N*C+1 <= i <= N*D) then (-1 + (i-N*C)*2/(N*(D-C)))
		else  1.0;

var u{i in 0..N} <= 1.0, >= -1.0, := v[i];
var w{0..N} := 0.0;

minimize f:
	-2*u[0]*u[1] + u[0]^2 + sum {i in 1..N-1} (-2*u[i]*u[i+1] + 2*u[i]^2) +
	u[N]^2 + sum {i in 0..N} 1/N^2*u[i]*w[i]+
	sum {i in 0..N} (-1/N^2*v[i]*u[i] - 2/N^2*v[i]*w[i]) + (v[1]-v[0])*u[0]
	+ sum {i in 1..N-1} (v[i-1]-2*v[i]+v[i+1])*u[i] + (v[N-1]-v[N])*u[N];
subject to cons1:
	0.5*u[0] + sum{i in 1..N-1} u[i] + 0.5*u[N] = 0.2*INT;
subject to cons2:
	u[0] - u[1] - 1/N^2*w[0] = 0;
subject to cons3{i in 1..N-1}:
	2*u[i] - u[i+1] - u[i-1] - 1/N^2*w[i] = 0;
subject to cons4:
	u[N] - u[N-1] - 1/N^2*w[N] = 0;

solve;
display f;
display u,w;	
