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

#   Source: an exercize for L. Watson course on LANCELOT in the Spring 1993.
#   B.Benhabib, R.G.Fenton and A.A.Goldberg, 
#   "Analytical trajectory optimization of seven degrees of freedom redundant
#   robot",  
#   Transactions of the Canadian Society for Mechanical Engineering,
#   vol.11(4), 1987, pp 197-200.

#   SIF input: Manish Sabu at Virginia Tech., Spring 1993.
#              Minor modifications by Ph. L. Toint, April 1993.

#   classification QOR2-MY-14-2

param XPOS := 4;
param YPOS := 4;
param HIGH := 2.356194;
param DOWN := -2.356194;
param THIN{1..7} := 0.0;

var TH{1..7};
var THI{1..7};

minimize f:
	sum {i in 1..7} (TH[i]-THI[i])^2;
subject to cons1:
	sum {i in 1..6} cos(TH[i])+0.5*cos(TH[7]) -XPOS = 0;
subject to cons2:
	sum {i in 1..6} sin(TH[i])+0.5*sin(TH[7]) -YPOS = 0;
subject to cons3{i in 1..7}:
	THI[i] = THIN[i];

option loqo_options "verbose=2 timing=1 iterlim=5000";
solve;
display f;
display TH, THI;
