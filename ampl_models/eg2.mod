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

#   classification OUR2-AN-1000-0

param N := 1000;

var x{1..N};

minimize f:
	sum {i in 1..N-1} sin(x[1] + x[i]^2 - 1.0) + 0.5*sin(x[N]^2);

solve;
display f;
display x;
