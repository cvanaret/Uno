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

#   classification OOR2-MY-6-8

param np := 4;
param x{1..np};
param y{1..np};

var v1 := -40.0;
var w1 := 5.0;
var d >= 1e-8, := 1.0;
var a >= 1.0, := 2.0;
var t >= 0.0, <= 6.2831852, := 1.5;
var r >= 0.39, := 0.75;

minimize f:
	(d+r)^2*acos(-( (a*d)^2 - (a*d+r)^2 +(d+r)^2)/(2*(d+r)*a*d))
	-(a*d+r)^2*acos(( (a*d)^2+ (a*d+r)^2 -(d+r)^2)/(2*(a*d+r)*a*d))
	+(d+r)*a*d*sin(acos(-( (a*d)^2 - (a*d+r)^2 +(d+r)^2)/(2*(d+r)*a*d)));
subject to cons1{i in 1..np}:
	(v1+a*d*cos(t)-x[i])^2 + (w1+a*d*sin(t)-y[i])^2 - (d+r)^2<= 0.0;
subject to cons2{i in 1..np}:
	(v1-x[i])^2 + (w1-y[i])^2 - (a*d+r)^2 >= 0.0;
data;
param: x y :=
1	1.0	0.0
2	0.0	1.0
3	0.0	-1.0
4	0.5	0.0
;

solve;
display f;
display v1, w1, a, d, r, t;
