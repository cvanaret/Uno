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
#   M.M. Makela,
#   "Nonsmooth optimization",
#   Ph.D. thesis, Jyvaskyla University, 1990

#   SIF input: Ph. Toint, Nov 1993.

#   classification  LLR2-AN-21-40

var x{1..20};
var u;

minimize f:
	u;
subject to cons1{i in 1..20}:
	-u+x[i] <= 0;
subject to cons2{i in 1..20}:
	-u-x[i] <= 0;

data;
var x:=
1	1
2	2
3	3
4	4
5	5
6	6
7	7
8	8
9	9
10	10
11	-11
12	-12
13	-13
14	-14
15	-15
16	-16
17	-17
18	-18
19	-19
20	-20;

solve;
display f;
display x;
