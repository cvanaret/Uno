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

#   classification QLR2-AN-V-V

param N:=1000;
param M:=N/2;
param NPLUS:=N/2;

var x{1..N} >= 0.1, <= 10.0, := 0.5;

minimize f:
	sum {i in 1..NPLUS} (x[i]+x[(2*i-1)-((2*i-1) div N)*N+1]
	+x[(3*i-1)-((3*i-1) div N)*N+1])^2*i/2 -
	sum {i in NPLUS+1..N} (x[i]+x[(2*i-1)-((2*i-1) div N)*N+1]
	+x[(3*i-1)-((3*i-1) div N)*N+1])^2*i/2 ;
subject to cons1{i in 1..M}:
	x[i]+2*x[(4*i-1)-((4*i-1) div N)*N+1] + 3*x[(5*i-1)-((5*i-1)
	div N)*N+1] = 6.0 ;

solve;
display f;
display x;


