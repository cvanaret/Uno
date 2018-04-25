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
#   J. Hald and K. Madsen,
#   "Combined LP and quasi-Newton methods for minimax optimization",
#   Mathematical Programming 20, pp. 49-62, 1981.

#   SIF input: Ph. Toint, Nov 1993.

#   classification  LOR2-AN-6-42

param y{i in 1..21} := -1.0+0.1*(i-1);
param ey{i in 1..21} := exp(y[i]);

var x{i in 1..5} := if (i==1) then 0.5 else 0.0;
var u := 0.0;

minimize f:
	u;
subject to cons1{i in 1..21}:
	(x[1]+y[i]*x[2])/(1.0+x[3]*y[i]+x[4]*y[i]^2+x[5]*y[i]^3)-u <= ey[i];
subject to cons2{i in 1..21}:
	-(x[1]+y[i]*x[2])/(1.0+x[3]*y[i]+x[4]*y[i]^2+x[5]*y[i]^3)-u <= -ey[i];

solve;
display f;
display x;
