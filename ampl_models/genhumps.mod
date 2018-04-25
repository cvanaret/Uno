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
#   Ph. Toint, private communication, 1997.

#   SDIF input: N. Gould and Ph. Toint, November 1997.

#   classification SUR2-AN-V-0

param zeta := 2;
param N:=5;

var x{i in 1..N} := if (i==1) then -506.0 else 506.2;

minimize f:
	sum {i in 1..N-1} ( sin (zeta*x[i])^2*sin(zeta*x[i+1])^2+0.05*(x[i]^2+x[i+1]^2) );

solve;
display f;
display x;
