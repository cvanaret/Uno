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

#  Source:
#  D. Tuyttens,
#  "A User's Guide to LSNNO, a Fortran subroutine for large scale
#  nonlinear optimization problems",
#  Report 89/00, Department of Mathemetics, FUNDP, Namur (Belgium), 1989.

#  SIF input: J-M COLLIN, July 1990.

#  classification ONR2-AY-5-4

param N:=5;
param x_init{1..5};
param x_lower{1..5} default -Infinity;
param x_upper{1..5} default Infinity;

var x{i in 1..N} >= x_lower[i], <= x_upper[i], := x_init[i];

minimize f:
	x[2]*exp(x[1]+x[3]) + x[3]^2*x[4]^2 + (x[3]-x[5])^2;
subject to cons1:
	x[1]+x[2]-10 = 0;
subject to cons2:
	-x[1]-x[3]+x[4] = 0;
subject to cons3:
	-x[2]+x[3]+x[5] = 0;
subject to cons4:
	-x[4]-x[5]+10 = 0;

data;
param x_init:=
1	4
2	6
3	2
4	6
5	4;
param:
	x_lower	x_upper:=
1	2	4
2	6	8
3	0	5;

solve;
display f;
display x;
