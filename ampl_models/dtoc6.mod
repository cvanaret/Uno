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

#   classification OOR2-AN-V-V

param n := 5001;

var x{1..n-1};
var y{1..n};

minimize f:
	sum {t in 1..n-1} ( 0.5*(exp(x[t])+y[t])^2 + 0.5*(x[t])^2 );
subject to cons1{t in 1..n-1}:
	y[t] - y[t+1] + exp(x[t]) = 0;
fix y[1] := 0.0;
solve;	
