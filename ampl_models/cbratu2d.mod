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

#   Source: problem 3 in
#   J.J. More',
#   "A collection of nonlinear model problems"
#   Proceedings of the AMS-SIAM Summer seminar on the Computational
#   Solution of Nonlinear Systems of Equations, Colorado, 1988.
#   Argonne National Laboratory MCS-P60-0289, 1989.

#   SIF input: Ph. Toint, Dec 1989.

#   classification NOR2-MN-V-V

param p := 23;
param lambda := 5.0;
param h := 1/(p-1);
param c := h^2/lambda;

var u{i in 1..p, j in 1..p} := 0.0;
var x{i in 1..p, j in 1..p} := 0.0;

minimize f:0;
subject to cons1{i in 2..p-1, j in 2..p-1}:
	(4*u[i,j]-u[i+1,j]-u[i-1,j]-u[i,j+1]-u[i,j-1]-c*exp(u[i,j])*cos(x[i,j])) = 0;
subject to cons2{i in 2..p-1, j in 2..p-1}:
	(4*x[i,j]-x[i+1,j]-x[i-1,j]-x[i,j+1]-x[i,j-1]-c*exp(u[i,j])*sin(x[i,j])) = 0;

fix {j in 1..p} u[1,j] := 0.0;
fix {j in 1..p} u[p,j] := 0.0;
fix {j in 1..p} x[1,j] := 0.0;
fix {j in 1..p} x[p,j] := 0.0;

fix {i in 2..p-1} u[i,p] := 0.0;
fix {i in 2..p-1} u[i,1] := 0.0;
fix {i in 2..p-1} x[i,p] := 0.0;
fix {i in 2..p-1} x[i,1] := 0.0;

solve;
display f;
display u,x;
	
