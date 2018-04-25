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

#   classification LOR2-AN-7-V

param m := 500;
param lower := -0.5;
param upper := 0.5;
param diff := upper-lower;
param h := diff/m;

var u  := 0.0;
var x{1..6} := 0.0;

minimize f:
	u;
subject to cons1{i in 0..m}:
	-(x[1]*exp(x[4]*(i*h+lower)) + x[2]*exp(x[5]*(i*h+lower)) + x[3]*exp(x[6]*(i*h+lower))) + u >= -1/(1+i*h+lower);
subject to cons2{i in 0..m}:
	(x[1]*exp(x[4]*(i*h+lower)) + x[2]*exp(x[5]*(i*h+lower)) + x[3]*exp(x[6]*(i*h+lower))) + u >= 1/(1+i*h+lower);

solve;
display f;
display x,u;

