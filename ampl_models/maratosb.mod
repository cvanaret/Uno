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
#   Ph. Toint.

#   SIF input: Ph. Toint, June 1990.

#   classification QUR2-AN-2-0

param invp := 0.000001;
param xinit{1..2};
var x{i in 1..2} := xinit[i];

minimize f:
	x[1]+(x[1]^2+x[2]^2-1)^2/invp;

data;
param xinit:= 1 0 2 0;

solve; display f; display x;
