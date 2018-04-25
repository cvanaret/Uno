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
#   M. Aljazzaf,
#   "Multiplier methods with partial elimination of constraints for
#   nonlinear programming",
#   PhD Thesis, North Carolina State University, Raleigh, 1990.

#   SDIF input: Ph. Toint, May 1990.

#   classification QQR2-AN-3-1

param N := 3;
param N1 := 2;
param Biga := 100.0;
param F := (Biga^2-1.0)/(N-1);
param F2 := (Biga^2-1.0)/(Biga*(N-1));
param A{i in 1..N} := Biga-(i-1)*F2;
param B{i in 1..N} := (i-1)*F+1.0;
var x{1..N} := 0.0, >= 0.0;

minimize f:
	A[1]*(x[1]-0.5)^2 + sum {i in 2..N1} A[i]*(x[i]+1.0)^2 + sum {i in N1+1..N} A[i]*(x[i]-1.0)^2;

subject to cons1:
	-B[1]*x[1]+B[1]+sum {i in 2..N1} B[i]*(x[i]-0.0)^2+sum {i in N1+1..N} B[i]*(x[i]-1.0)^2 = 0;

display A, B, F, F2, f;
solve;
display f;
display x;
