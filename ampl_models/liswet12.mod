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
#   W. Li and J. Swetits,
#   "A Newton method for convex regression, data smoothing and
#   quadratic programming with bounded constraints",
#   SIAM J. Optimization 3 (3) pp 466-488, 1993.

#   SIF input: Nick Gould, August 1994.

#   classification QLR2-AN-V-V

param N:=10000;
param K:=2;
param B{i in 0..K} := if (i=0) then 1 else B[i-1]*i;
param C{i in 0..K} := if (i=0) then 1 else (-1)^i*B[K]/(B[i]*B[K-i]);
param T{i in 1..N+K} := (i-1)/(N+K-1);
param pi:=3.1415;
var x{1..N+K} := 0.0;

minimize f:
	sum {i in 1..N+K} -(cos(4*pi*T[i])+0.1*sin(i))*x[i] -
	sum {i in 1..N+K} 0.5*(cos(4*pi*T[i])+0.1*sin(i))^2+
	sum {i in 1..N+K} 0.5*x[i]^2;
subject to cons1{j in 1..N}:
	sum {i in 0..K} C[i]*x[j+K-i] >= 0;

solve;
display f;
display x;
