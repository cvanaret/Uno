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

#   Source: Problem 8, eqs (8.6)--(8.9) in
#   J.J. More',
#   "A collection of nonlinear model problems"
#   Proceedings of the AMS-SIAM Summer seminar on the Computational
#   Solution of Nonlinear Systems of Equations, Colorado, 1988.
#   Argonne National Laboratory MCS-P60-0289, 1989.

#   SIF input: Ph. Toint, Dec 1989.
#              minor correction by Ph. Shott, Jan 1995 and F Ruediger, Mar 1997.

#   classification NOR2-MN-V-V

param n :=2500;
param pem := 1.0;
param peh := 5.0;
param d := 0.135;
param b := 0.5;
param beta := 2.0;
param gamma := 25.0;
param h := 1/(n-1);
param cu1 := -h*pem;
param cui1 := 1/(h^2*pem)+1/h;
param cui := -1/h - 2/(h^2*pem);
param ct1 := -h*peh;
param cti1 := 1/(h^2*peh)+1/h;
param cti := -beta -1/h - 2/(h^2*peh);

var t{1..n} >= 0.0000001, := 1.0;
var u{1..n} >= 0.0, := 1.0;

minimize f: 0;
subject to cons1:
	(cu1*u[2]-u[1]+h*pem) = 0;
subject to cons2:
	(ct1*t[2]-t[1]+h*peh) = 0;
subject to cons3{i in 2..n-1}: 
	(-d*u[i]*exp(gamma-gamma/t[i])+(cui1)*u[i-1] + cui*u[i] + u[i+1]/(h^2*pem)) = 0;
subject to cons4{i in 2..n-1}:
	(b*d*u[i]*exp(gamma-gamma/t[i])+(cti1)*t[i-1] + cti*t[i] + t[i+1]/(h^2*peh)) = 0;
subject to cons5:
	(u[n]-u[n-1]) = 0;
subject to cons6:
	(t[n]-t[n-1]) = 0;

solve;
display f;
display t, u;
