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

#   Source: Nick Gould

#   SIF input: Nick Gould, September 1997.

#   classification SUR2-AN-V-0

param N:=10000;
param K:=10;
param sc := 12.0;
param scale{i in 1..N} := exp((i-1)*sc/(N-1));
var x{i in 1..N} := 0.0001*scale[i]/(N+1);
var y{i in 1..N-K} = sum {j in i..i+K} x[j]*scale[j];
var z{i in N-K+1..N} = sum {j in i..N} x[j]*scale[j];

minimize f:
	sum {i in 1..N-K} (y[i]*(y[i]*(y[i]^2-20)-0.1)) + sum {i in N-K+1..N} (z[i]*(z[i]*(z[i]^2-20)-0.1));

solve; display f; display x;
