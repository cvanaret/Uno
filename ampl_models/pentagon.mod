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
#   M.J.D. Powell,
#   " TOLMIN: a Fortran package for linearly constrained
#   optimization problems",
#   Report DAMTP 1989/NA2, University of Cambridge, UK, 1989.

#   SIF input: Ph. Toint, May 1990.

#   classification OLR2-AY-6-15

param pi := 4*atan(1);

param s{j in 0..4} := sin(2*pi*j/5);
param c{j in 0..4} := cos(2*pi*j/5);

var x{1..3};
var y{1..3};

minimize f:
	1/( (x[1]-x[2])^2+(y[1]-y[2])^2 )^8 +
	1/( (x[1]-x[3])^2+(y[1]-y[3])^2 )^8 +
	1/( (x[3]-x[2])^2+(y[3]-y[2])^2 )^8;
subject to cons{i in 1..3, j in 0..4}:
	c[j]*x[i] + s[j]*y[i] <= 1.0;

data;
var x:= 
1 -1
2 0
3 1;

var y := 
1 0
2 -1
3 1;

solve;
display f;
display x,y;
