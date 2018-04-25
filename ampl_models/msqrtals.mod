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

#   Source:  problem 201 (p. 93) in
#   A.R. Buckley,
#   "Test functions for unconstrained minimization",
#   TR 1989CS-3, Mathematics, statistics and computing centre,
#   Dalhousie University, Halifax (CDN), 1989.

#   SIF input: Ph. Toint, Dec 1989.

#   classification SUR2-AN-V-V

param P := 32;
param N := P^2;
param B{i in 1..P, j in 1..P} := sin ( ((i-1)*P + j)^2 );
param A{i in 1..P, j in 1..P} := sum {k in 1..P} B[i,k]*B[k,j];

var x{i in 1..P, j in 1..P} := 0.2*B[i,j];

minimize f:
	sum {i in 1..P, j in 1..P} (sum {t in 1..P} x[i,t]*x[t,j]
	-A[i,j])^2;

solve;
display f;
display x;
