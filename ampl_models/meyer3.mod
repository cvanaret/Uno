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

#   Source:  Problem 28 in
#   J.J. More', B.S. Garbow and K.E. Hillstrom,
#   "Testing Unconstrained Optimization Software",
#   ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.

#   See also Buckley #29 (p. 73).

#   SIF input: Ph. Toint, Dec 1989.

#   classification SUR2-RN-3-0

param N:=3;
param M:=16;

param x_init{1..N};
var x{i in 1..N} := x_init[i];

param t{i in 1..M} := 45+5*i;
param y{1..M};

minimize f:
	sum {i in 1..M} (x[1]*exp(x[2]/(t[i]+x[3]))-y[i])^2;

data;
param x_init:= 1 0.02 2 4000 3 250;

param y:=
1	34780
2	28610
3	23650
4	19630
5	16370
6	13720
7	11540
8	9744
9	8261
10	7030
11	6005
12	5147
13	4427
14	3820
15	3307
16	2872;

solve;
display f;
display x;
