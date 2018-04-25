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

#   Source:  problem 8, eqs (8.10)--(8.11) in
#   J.J. More',
#   "A collection of nonlinear model problems"
#   Proceedings of the AMS-SIAM Summer seminar on the Computational
#   Solution of Nonlinear Systems of Equations, Colorado, 1988.
#   Argonne National Laboratory MCS-P60-0289, 1989.

#   SIF input: Ph. Toint, Dec 1989.

#   classification NOR2-MN-V-V

param n := 1000;

param pe := 5.0;
param d := 0.135;
param b := 0.5;
param gamma := 25.0;

param h := 1/(n-1);
param ct1 := -h*pe;
param cti1 := 1/h + 1/(h^2*pe);
param cti := -1/h-2/(h^2*pe);

var t{1..n} >= 0.0000001, := 1.0;

minimize f: 0;
subject to cons1:
	(ct1*t[2]-t[1]+h*pe) = 0;
subject to cons2{i in 2..n-1}:
	(d*(b+1-t[i])*exp(gamma-gamma/t[i])+cti1*t[i-1]+cti*t[i]+t[i+1]/(h^2*pe)) = 0;
subject to cons3:
	(t[n]-t[n-1]) = 0;

solve;
display f;
display t;
