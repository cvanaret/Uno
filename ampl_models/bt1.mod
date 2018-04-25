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

#   Source: problem 1 in
#   P.T. Boggs and J.W. Tolle,
#   "A strategy for global convergence in a sequential 
#    quadratic programming algorithm",
#   SINUM 26(3), pp. 600-623, 1989.
 
#   SIF input: Ph. Toint, June 1993.

#   classification QQR2-AN-2-1

param xinit{1..2};
var x{1..2};
minimize f:
	100*x[1]^2+100*x[2]^2-x[1]-100;
subject to cons1:
	x[1]^2+x[2]^2-1.0=0;

data;
param xinit:= 1 0.08 2 0.06;

solve;
display f;
display x;


