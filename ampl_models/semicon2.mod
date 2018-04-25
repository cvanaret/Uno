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

#   Source: problem 10 in
#   J.J. More',
#   "A collection of nonlinear model problems"
#   Proceedings of the AMS-SIAM Summer seminar on the Computational
#   Solution of Nonlinear Systems of Equations, Colorado, 1988.
#   Argonne National Laboratory MCS-P60-0289, 1989.

#   SIF input: Ph. Toint, Dec 1989.

#   classification NOR2-AN-V-V

param A := -0.00009;
param B := 0.00001;

param N := 1000;
param LN := 900;

param UA := 0.0;
param UB := 700.0;

param CA := 1.0D12;
param CB := 1.0D13;
param BETA := 40.0;

param LAMBDA := 0.2;

param H := (B-A)/(N+1);
param LB := LAMBDA*BETA;
param LUA := LAMBDA*UA;
param LUB := LAMBDA*UB;
param ULW := LUA-5;
param UUP := LUB+5;

var u{i in 0..N+1} <= UUP, >= ULW, := if (i == 0) then LUA
			else if (i==N+1) then LUB else 0.0;

minimize f: 0;
subject to cons1{i in 1..LN}:
	(u[i-1]-2*u[i]+u[i+1]-LAMBDA*H^2*CA + LAMBDA*H^2*CA*exp(-LB*(u[i]-LUA))) = 0;
subject to cons2{i in LN+1..N}:
	(u[i-1]-2*u[i]+u[i+1]+LAMBDA*H^2*CB - LAMBDA*H^2*CB*exp(LB*(u[i]-LUB))) = 0;
fix u[0] := LUA;
fix u[N+1] := LUB;

solve;
display f;
display u;
