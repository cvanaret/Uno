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

#   Source: problem 11 in
#   J.J. More', B.S. Garbow and K.E. Hillstrom,
#   "Testing Unconstrained Optimization Software",
#   ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.

#   See also Buckley#27
#   SIF input: Ph. Toint, Dec 1989.

#   classification SUR2-MN-3-0

param N:=3;
param M:=99;

param t{i in 1..M} := i/100;
param y{i in 1..M} := 25+ (-50*log(t[i]))^(2/3);

param x_init{1..N};
var x{i in 1..N} := x_init[i];

minimize f:
	sum {i in 1..M} (exp(abs(y[i]-x[2])^x[3]/(-x[1]))-t[i])^2;

data;
param x_init := 1 5 2 2.5 3 0.15;

solve;
display f;
display x;
