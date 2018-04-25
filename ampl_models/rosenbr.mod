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

#   Source:  problem 1 in
#   J.J. More', B.S. Garbow and K.E. Hillstrom,
#   "Testing Unconstrained Optimization Software",
#   ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.

#   SIF input: Ph. Toint, Dec 1989, Corrected July 1993.

#   classification SUR2-AN-2-0

param xinit{1..2};
var x{i in 1..2} := xinit[i];

minimize f:
	(x[2]-x[1]^2)^2/0.01+(x[1]-1)^2;

data;
param xinit:= 1 -1.2 2 1.0;

solve; display f; display x;
