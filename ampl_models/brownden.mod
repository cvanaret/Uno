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

#   Source: Problem 16 in
#   J.J. More', B.S. Garbow and K.E. Hillstrom,
#   "Testing Unconstrained Optimization Software",
#   ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.

#   See also Buckley#30
#   SIF input: Ph. Toint, Dec 1989.

#   classification SUR2-AN-4-0

param M:=20;
param N:=4;

param x_init{1..4};
var x{i in 1..4} := x_init[i];

param t{i in 1..M} := i/5;

minimize f:
	sum {i in 1..M} ( (x[1]+t[i]*x[2]-exp(t[i]))^2 + (x[3]+x[4]*sin(t[i])-cos(t[i]))^2 )^2;

data;
param x_init:= 1 25 2 5 3 -5 4 -1;

solve;
display f;
display x;
