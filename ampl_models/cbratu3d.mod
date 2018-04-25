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

#   Source: Problem 3 in
#   J.J. More',
#   "A collection of nonlinear model problems"
#   Proceedings of the AMS-SIAM Summer seminar on the Computational
#   Solution of Nonlinear Systems of Equations, Colorado, 1988.
#   Argonne National Laboratory MCS-P60-0289, 1989.

#   SIF input: Ph. Toint, Dec 1989.

#   classification NOR2-MN-V-V

param p := 10;
param lambda := 6.80812;
param h := 1/(p-1);
param c := h^2/lambda;

var u{i in 1..p, j in 1..p, k in 1..p} := 0.0;
var x{i in 1..p, j in 1..p, k in 1..p} := 0.0;

minimize f:0;
subject to cons1{i in 2..p-1, j in 2..p-1, k in 2..p-1}:
	(6*u[i,j,k]-u[i+1,j,k]-u[i-1,j,k]-u[i,j+1,k]-u[i,j-1,k]-u[i,j,k-1]-u[i,j,k+1]-c*exp(u[i,j,k])*cos(x[i,j,k])) = 0;
subject to cons2{i in 2..p-1, j in 2..p-1, k in 2..p-1}:
	(6*x[i,j,k]-x[i+1,j,k]-x[i-1,j,k]-x[i,j+1,k]-x[i,j-1,k]-u[i,j,k-1]-u[i,j,k+1]-c*exp(u[i,j,k])*sin(x[i,j,k])) = 0;

fix {j in 1..p, k in 1..p} u[1,j,k] := 0.0;
fix {j in 1..p, k in 1..p} u[p,j,k] := 0.0;
fix {j in 1..p, k in 1..p} x[1,j,k] := 0.0;
fix {j in 1..p, k in 1..p} x[p,j,k] := 0.0;

fix {i in 2..p-1, k in 1..p} u[i,p,k] := 0.0;
fix {i in 2..p-1, k in 1..p} u[i,1,k] := 0.0;
fix {i in 2..p-1, k in 1..p} x[i,p,k] := 0.0;
fix {i in 2..p-1, k in 1..p} x[i,1,k] := 0.0;

fix {i in 2..p-1, j in 2..p-1} u[i,j,1] := 0.0;
fix {i in 2..p-1, j in 2..p-1} u[i,j,p] := 0.0;
fix {i in 2..p-1, j in 2..p-1} x[i,j,1] := 0.0;
fix {i in 2..p-1, j in 2..p-1} x[i,j,p] := 0.0;

solve;
display f;
display u,x;
	
