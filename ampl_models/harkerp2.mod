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
#   P. T. Harker and J.-S. Pang,
#   "A damped Newton method for the linear complementarity problem",
#   in 'Allgower and Georg: Computational solution of nonlinear
#   systems of equations', AMS lectures in Applied Mathematics 26,
#   AMS, Providence, Rhode Island, USA, pp 265-284.

#   SIF input: Nick Gould, July 1993.
#   classification QBR2-AN-V-V

param N:=100;
var x{i in 1..N} := i, >= 0.0;
minimize f:
	sum {i in 1..N} -1*x[i]^2*0.5 +
	sum {i in 1..N} -x[i] +
	( sum {i in 1..N} x[i] )^2 +
	sum {j in 2..N} 2*(sum {i in j..N} x[i])^2;

solve;
