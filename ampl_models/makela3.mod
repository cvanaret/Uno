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
#   M.M. Makela,
#   "Nonsmooth optimization",
#   Ph.D. thesis, Jyvaskyla University, 1990

#   SIF input: Ph. Toint, Nov 1993.

#   classification  LQR2-AN-21-20

var x{i in 1..20} := if (i <= 10) then i else -i;
var u;

minimize f:
	u;
subject to cons{i in 1..20}:
	-u+x[i]^2 <= 0;
solve;
