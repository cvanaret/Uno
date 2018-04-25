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
#   an example in a talk by W.K. Zhang and C. Fleury, LLN, 1994.

#   SIF input: Ph. Toint, November 1994

#   classification LOR2-MN-5-1

param num{1..5};
var x{1..5} >= 0.000001, := 1.0;

minimize f:
	0.0624*sum{i in 1..5} x[i];
subject to cons1:
	sum {i in 1..5} num[i]/x[i]^3 - 1 <= 0;

data;
param num:=
1	61.0
2	37.0
3	19.0
4	7.0
5	1.0;

solve;
display f;
display x;


