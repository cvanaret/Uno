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

#   Source:  problem 4.2 in
#   Y. Xiao and F. Zhou,
#   "Non-monotone trust region methods with curvilinear path
#   in unconstrained optimization",
#   Computing, vol. 48, pp. 303-317, 1992.

#   SIF input: F Facchinei, M. Roma and Ph. Toint, June 1994

#   classification SUR2-AN-2-0

param c:=10D-4;
param xinit{1..2};
var x{i in 1..2} := xinit[i];
minimize f:
	(x[2]-sin(x[1]))^2/c + x[1]^2/4;

data;
param xinit:= 1 4.712389 2 -1.0;

solve; display f; display x;
