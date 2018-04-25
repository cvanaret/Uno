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
#   R. Byrd,
#   Private communication, Chicago, 1992.
 
#   SIF input: Ph. Toint, November 1992.
 
#   classification LQR2-AN-3-2

param xinit{1..3};
var x{i in 1..3} := xinit[i];

minimize f:
	-x[1]-x[2]-x[3];
subject to cons1:
	-9.0+x[1]^2+x[2]^2+x[3]^2 = 0;
subject to cons2:
	-9.0+(x[1]-1.0)^2+x[2]^2+x[3]^2 = 0;

data;
param xinit:=
1	5.0
2	0.0001
3	-0.0001;

solve;
display f;
display x;
