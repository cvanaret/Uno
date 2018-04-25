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

#   Source: K. Oettershagen
#   "Ein superlinear konvergenter algorithmus zur losung 
#    semi-infiniter optimierungsproblem",
#    Ph.D thesis, Bonn University, 1982

#   SIF input: Nick Gould, February, 1994.

#   classification LLR2-AN-4-V

param M:=500;
param lower := 0.0;
param upper := 1.0;
param diff := upper-lower;
param h:=diff/M;

var u;
var x{1..3};

minimize f:
	u;
subject to cons1{i in 0..M}:
	u-sin(i*h+lower)-x[1]-(i*h+lower)*x[2] - (i*h+lower)^2*x[3] >= 0;
subject to cons2{i in 0..M}:
	u+sin(i*h+lower)+x[1]+(i*h+lower)*x[2] + (i*h+lower)^2*x[3] >= 0;
	
solve; display f; display x;
