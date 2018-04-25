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
#   A Griewank and Ph. Toint,
#   "Partitioned variable metric updates for large structured
#   optimization problems",
#   Numerische Mathematik 39:429-448, 1982.

#   SIF input: Ph. Toint, Dec 1989.

#   classification OXR2-MY-V-0

param P:=125;
param h00 := 1;
param slopej := 4;
param slopei := 8;

param ston:= slopei/(P-1);
param wtoe:= slopej/(P-1);
param h01:= h00+slopej;
param h10:= h00+slopei;

var x{i in 1..P,j in 1..P} := 0;

minimize f:
	sum {i in 1..P-1, j in 1..P-1}
	sqrt((P-1)^2*((x[i,j]-x[i+1,j+1])^2+(x[i+1,j]-x[i,j+1])^2)/2+1)/(P-1)^2;
subject to cons1{j in 1..P}:
	x[1,j] = (j-1)*wtoe+h00;
subject to cons2{j in 1..P}:
	x[P,j] = (j-1)*wtoe+h10;
subject to cons3{i in 2..P-1}:
	x[i,P] = (i-1)*ston+h01;
subject to cons4{i in 2..P-1}:
	x[i,1] = (i-1)*ston+h00;

solve;
display f;
display x;
