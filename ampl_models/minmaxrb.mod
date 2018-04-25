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
#   J. Hald and K. Madsen
#   "Combined LP and quasi-Newton methods for minmax optimization"
#   Mathematical Programming, vol.20, p. 49-62, 1981.

#   SIF input: Ph. Toint, April 1992.

#   classification LQR2-AN-3-4

param xinit{1..2};
var x{i in 1..2} := xinit[i];
var u := 1.0;
minimize f:
	u;
subject to cons1:
	u-10*x[2]+10*x[1]^2 >= 0;
subject to cons2:
	u+10*x[2]-10*x[1]^2 >= 0;
subject to cons3:
	u+x[1]-1 >= 0;
subject to cons4:
	u-x[1]+1 >= 0;

data;
param xinit:= 1 -1.2 2 1.0;

solve;
display f;
display x;
display u;
