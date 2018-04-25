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
#   "The OPTIMA user manual (issue No.8, p. 49)",
#   Numerical Optimization Centre, Hatfield Polytechnic (UK), 1989.

#   SIF input: Ph. Toint, May 1990.

#   classification NOR2-AY-25-25

param N:=25;
var x{1..N} := 1.0;

minimize f: 0;
subject to cons{i in 2..N-1}: 
	( x[i]*(x[i-1]-x[i+1])+x[i]-x[13]+1 ) = 0;
subject to cons2: 
	(x[1]-x[13]+1-x[1]*x[2]) = 0;
subject to cons3:
	(x[N]-x[13]+1+x[N-1]*x[N]) = 0;

solve; display f; display x;
