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
#   Ph. Toint

#   SIF input: Ph. Toint, Dec 1989.

#   classification SUR2-AN-V-0

param P := 3;
param N := P^2;
param B{i in 1..P, j in 1..P} := if (i==3 && j==1) then 0 else sin (
				 ((i-1)*P + j)^2 );
param A{i in 1..P, j in 1..P} := sum {k in 1..P} B[i,k]*B[k,j];

var x{i in 1..P, j in 1..P} := if (i==3 && j==1) then -0.8*sin(((i-1)*P + j)^2 ) else 0.2*B[i,j];

minimize f:
	sum {i in 1..P, j in 1..P} (sum {t in 1..P} x[i,t]*x[i,j]
	-A[i,j])^2;

option loqo_options "verbose=2 timing=1 iterlim=400 sigfig=5 inftol=0.00001";
solve;
display f;
display x;
