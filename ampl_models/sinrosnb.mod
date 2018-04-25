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

#   Source: a modification of an original idea by
#   Ali Bouriacha, private communication.

#   SIF input: Nick Gould and Ph. Toint, October, 1993.

#   classification OQR2-AN-V-V

param N:=1000;
param pi:=4*atan(1);

var x{i in 1..N} := if (i==1) then 1 else 1;

minimize f:
	sin(x[1]-1+1.5*pi)
	+ sum {i in 2..N} 100*sin(-x[i]+1.5*pi+x[i-1]^2)
;
subject to cons1{i in 2..N}:
	2*pi >= pi-x[i]+x[i-1]^2 >= 0;
subject to cons2:
	-pi <= x[1] <= pi;
;

solve;
display f;
display x;
