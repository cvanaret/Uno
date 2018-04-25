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

#   classification QBR2-AN-V-0

param N:=10000;
param M:=N/2;
param Nplus :=N/2;
var x{1..N} <= 10.0, >= 0.1, := 0.5;

minimize f:
	sum {i in 1..Nplus} 0.5*i*(x[i]+x[ ((2*i-1) mod N) + 1] + x[ ((3*i-1) mod N) +1 ] )^2 -
	sum {i in Nplus+1..N} 0.5*i*(x[i]+x[ ((2*i-1) mod N) + 1] + x[ ((3*i-1) mod N) +1 ] )^2 ;

solve; display f; display x;
