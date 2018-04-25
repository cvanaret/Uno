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

#   classification OUR2-AN-V-0

param n := 1000;
param alpha := 0.5;

var x{i in 1..n} := i/(n+1);

minimize f:
	sum {i in 1..n} (x[i]) +
	sum {i in 2..n-1} alpha*cos(2*x[i]-x[n]-x[1]);

solve;
display f;
display x;
