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
#   K.M. Irani, M.P. Kamat, C.J. Ribbens, H.F.Walker and L.T. Watson,
#   "Experiments with conjugate gradient algoritms for homotopy curve
#    tracking" ,
#   SIAM Journal on Optimization, May 1991, pp. 222-251, 1991.

#   SIF input: Ph. Toint, May 1990.

#   classification NOR2-AN-V-V

param N := 5000;

var x{0..N+1} := 1.0;

minimize f: 0;
subject to cons{i in 1..N}:
	(-0.05*(x[i] + x[i+1] + x[i-1]) + atan( sin( (i mod 100)*x[i] ) )) = 0;

fix x[0] := 0.0;
fix x[N+1] := 0.0;

solve; display x;

