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

#   Source: setting the boundary free in 
#   A Griewank and Ph. Toint,
#   "Partitioned variable metric updates for large structured
#   optimization problems",
#   Numerische Mathematik 39:429-448, 1982.

#   SIF input: Ph. Toint, November 1991.

#   classification OUR2-MY-V-0

param p := 32;

param h00 := 1.0;
param slopej := 4.0;
param slopei := 8.0;

param scale := (p-1)^2;

param ston := slopei/(p-1);
param wtoe := slopej/(p-1);
param h01 := h00+slopej;
param h10 := h00+slopei;

var x{1..p,1..p};

minimize f:
	sum {i in 1..p-1, j in 1..p-1} sqrt(0.5*(p-1)^2*((x[i,j]-x[i+1,j+1])^2+(x[i+1,j]-x[i,j+1])^2)+1.0)/scale +
	(sum {j in 1..p, i in 1..p} x[i,j])^2/p^4;

let {j in 1..p} x[1,j] := (j-1)*wtoe+h00;
let {j in 1..p} x[p,j] := (j-1)*wtoe+h10;
let {i in 2..p-1} x[i,p] := (i-1)*ston+h00;
let {i in 2..p-1} x[i,1] := (i-1)*ston+h01;
let {i in 2..p-1,j in 2..p-1} x[i,j] := 0.0;

solve;
display f;
display x;
