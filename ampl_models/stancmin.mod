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
#   I.M. Stancu-Minasian
#   Private communication, 1992.

#   SDIF input: Ph. Toint, October 1992.

#   classification OLI2-AY-3-2

var x{1..3} := 50.0, >= 0;

minimize f:
	-(6*x[2]+3*x[1]+2*x[3]-11)/(x[1]+4*x[2]+x[3]+1);
subject to cons1:
	3*x[1]+4*x[2]+x[3]-2 <= 0;
subject to cons2:
	x[1]+4*x[2]+x[3]-1 <= 0; 

solve; display f; display x;
