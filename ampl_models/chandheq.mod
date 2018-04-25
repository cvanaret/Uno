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

#   Source: problem 4 in
#   J.J. More',
#   "A collection of nonlinear model problems"
#   Proceedings of the AMS-SIAM Summer seminar on the Computational
#   Solution of Nonlinear Systems of Equations, Colorado, 1988.
#   Argonne National Laboratory MCS-P60-0289, 1989.

#   SIF input: Ph. Toint, Dec 1989.

#   classification NOR2-RN-V-V

param n := 100;
param c := 1;
param x{i in 1..n} := i/n;
param w{i in 1..n} := 1/n;

var h{i in 1..n} := 1.0, >= 0.0;

minimize f: 0;
subject to cons{i in 1..n}:
	sum {j in 1..n} -0.5*c*w[j]*x[i]/(x[i]+x[j])*h[i]*h[j] + h[i] = 1.0;

solve;
display f;
display h;
