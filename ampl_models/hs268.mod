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
#   K. Schittkowski
#   "More Test Examples for Nonlinear Programming Codes"
#   Springer Verlag, Berlin, Lecture notes in economics and 
#   mathematical systems, volume 282, 1987

#   SIF input: Michel Bierlaire and Annick Sartenaer, October 1992.
#              minor correction by Ph. Shott, Jan 1995.

#   classification QLR2-AN-5-5

param D{1..5, 1..5};
param B{1..5};
var x{1..5} := 1.0;

minimize f:
	14463.0 + sum {i in 1..5, j in 1..5} D[i,j]*x[i]*x[j]
	+ -2*sum {i in 1..5} (B[i]*x[i]);
subject to cons1:
	-sum {i in 1..5} x[i] + 5>= 0;
subject to cons2:
	10*x[1]+10*x[2]-3*x[3]+5*x[4]+4*x[5] -20 >= 0;
subject to cons3:
	-8*x[1]+x[2]-2*x[3]-5*x[4]+3*x[5] + 40 >= 0;
subject to cons4:
	8*x[1]-x[2]+2*x[3]+5*x[4]-3*x[5] -11>= 0;
subject to cons5:
	-4*x[1]-2*x[2]+3*x[3]-5*x[4]+x[5] +30>= 0;

data;
param B:=
1	-9170
2	17099
3	-2271
4	-4336
5	-43;
param D:
	1	2	3	4	5:=
1	10197	-12454	-1013	1948	329
2	-12454	20909	-1733	-4914	-186
3	-1013	-1733	1755	1089	-174
4	1948	-4914	1089	1515	-22
5	329	-186	-174	-22	27;

solve; display f; display x;
