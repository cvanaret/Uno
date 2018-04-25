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
#   "Ein superlinear knonvergenter algorithmus zur losung 
#    semi-infiniter optimierungsproblem",
#    Ph.D thesis, Bonn University, 1982

#   SIF input: Nick Gould, February, 1994.

#   classification LOR2-AN-3-V

param M := 500;
param LOWER := -0.5;
param UPPER := 0.5;
param DIFF := UPPER-LOWER;
param H := DIFF/M;

var u := 0;
var x1 := 0;
var x2 := 0;

minimize f:
	u;
subject to cons1{i in 0..M}:
	u + 1/(1+i*H+LOWER) - x1*exp(x2*(i*H+LOWER))>= 0;
subject to cons2{i in 0..M}:
	u - 1/(1+i*H+LOWER) + x1*exp(x2*(i*H+LOWER))>= 0;

solve;
display f;
display u,x1,x2;

