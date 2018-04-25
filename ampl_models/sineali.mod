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

#   Source: on original idea by
#   Ali Bouriacha, private communication.

#   SIF input: Nick Gould and Ph. Toint, October, 1993.

#   classification OBR2-AN-V-0

param N:=20;
param pi := 3.1415926535;
var x{i in 1..N} >= if (i == 1) then -1.5*pi else sqrt(pi)-2*pi,
		 <= if (i == 1) then 0.5*pi else sqrt(pi),
		 := 0.0;

minimize f:
	sin(x[1]-1)
	+ sum {i in 2..N} 100*sin(x[i]-x[i-1]^2)
	;
option loqo_options "verbose=2 timing=1 inftol=0.000001 sigfig=6";
solve;
display f;
display x;
