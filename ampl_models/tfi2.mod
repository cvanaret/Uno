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
#   Y. Tanaka, M. Fukushima, T. Ibaraki,
#   "A comparative study of several semi-infinite nonlinear programming
#   algorithms",
#   EJOR, vol. 36, pp. 92-100, 1988.

#   SIF input: Ph. Toint, April 1992.

#   classification LLR2-AN-3-V

param M := 10000;
param h := 1/M;
var x{1..3};

minimize f:
	x[1]+0.5*x[2]+x[3]/3;
subject to cons1{i in 0..M}:
	-x[1]-i*h*x[2]-(i*h)^2*x[3]+tan(i*h) <= 0;

solve; display f; display x;
