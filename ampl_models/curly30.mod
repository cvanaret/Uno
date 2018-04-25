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
param K:=30;

var x{1..N} := 0.0001/(N+1);
var Q{i in 1..N} = if (i <= N-K) then (sum {j in i..i+K} x[j]) else (sum {j in i..N} x[j]);

minimize f:
	sum {i in 1..N} Q[i]*(Q[i]*(Q[i]^2-20)-0.1);

solve; display f; display x;
