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
#   A. R. Conn, N. I. M. Gould and Ph. L. Toint,
#   "LANCELOT: a Fortran package for large-scale nonlinear optimization
#    (Release A)", Springer Series in Computational Mathematics 17,
#   Springer Verlag, 1992

#   SIF input: A. R. Conn, Nick Gould and Ph. L. Toint, June 1994.

#   classification OOR2-AY-V-V

param n := 100;

var y;
var x{i in 1..n} >= -1.0, <= i, := 0.5;

minimize f:
	0.5*((x[1]-x[n])*x[2] + y)^2;
subject to consq{i in 1..n-1}:
	y + x[1]*x[i+1] + (1+2/i)*x[i]*x[n] <= 0.0;
subject to conss{i in 1..n}:
	0.5 >= (sin(x[i]))^2; # >= 0.0;
subject to eq:
	(x[1]+x[n])^2 = 1.0;

solve;
