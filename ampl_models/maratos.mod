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
#   A.A. Brown and M. Bartholomew-Biggs,
#   "Some effective methods for unconstrained optimization based on
#   the solution of ordinary differential equations",
#   Technical Report 178, Numerical Optimization Centre, Hatfield
#   Polytechnic, (Hatfield, UK), 1987.

#   SIF input: Nick Gould, June 1990.

#   classification QQR2-AN-2-1

param tau := 0.000001;
param xinit{1..2};
var x{i in 1..2} := xinit[i];

minimize f:
	-x[1]-tau+tau*(x[1]^2+x[2]^2); 
subject to cons1:
	-1.0+x[1]^2+x[2]^2 = 0;

data;
param xinit := 1 1.1 2 0.1;

solve;
display f;
display x;

