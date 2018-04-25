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

param P:=72;
param lambda := 6.80812;
param h:= 1/(P-1);
param c:= h^2*lambda;

var u{1..P,1..P};

minimize f: 0;
subject to cons{i in 2..P-1, j in 2..P-1}:
	(4*u[i,j]-u[i+1,j]-u[i-1,j]-u[i,j+1]-u[i,j-1]-c*exp(u[i,j])) = 0;
fix {j in 1..P} u[1,j] := 0;
fix {j in 1..P} u[P,j] := 0;
fix {i in 2..P-1} u[i,P] := 0;
fix {i in 2..P-1} u[i,1] := 0;

solve;
display f;
display u;
