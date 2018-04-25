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
#   "LANCELOT, A Fortran Package for Large-Scale Nonlinear Optimization
#   (Release A)"
#   Springer Verlag, 1992.

#   SIF input: N. Gould and Ph. Toint, June 1994.

#   classification OBR2-AY-3-0

	var x1;
	var x2 >= -1.0 ,  <= 1.0;
	var x3 >= 1.0 ,  <= 2.0;

minimize obj:
	x1^2 + (x2*x3)^4 + x1*x3 + x2*sin(x1+x3) + x2;

solve;
	display x1;
	display x2;
	display x3;
display obj;
