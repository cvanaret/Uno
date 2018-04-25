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

#   Source: Problem 12 in
#   J.J. More', B.S. Garbow and K.E. Hillstrom,
#   "Testing Unconstrained Optimization Software",
#   ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.

#   See also Buckley#BOX663
#   SIF input: Ph. Toint, Dec 1989.

#   classification SUR2-AN-3-0

param M:=10;

param x_init{1..3};
var x{i in 1..3} := x_init[i];

param t{i in 1..M} := 0.1*i;

minimize f:
	sum {i in 1..M} (exp(-t[i]*x[1])-exp(-t[i]*x[2])-x[3]*exp(-t[i])+x[3]*exp(-10*t[i]))^2;

data;
param x_init := 1 0 2 10 3 1.0;

solve;
display f;
display x;
