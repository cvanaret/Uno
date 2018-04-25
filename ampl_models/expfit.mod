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
#   A.R. Conn, N. Gould and Ph.L. Toint,
#   "LANCELOT, a Fortran package for large-scale nonlinear optimization",
#   Springer Verlag, FUNDP, 1992.

#   SIF input: Ph. Toint, Jan 1991.

#   classification SUR2-AN-2-0

param p:=10;
param h:=0.25; 
var alpha;
var beta;
minimize f:
	sum {i in 1..p} (alpha*exp(i*h*beta)-i*h)^2;

solve;
display f;
display alpha, beta;

