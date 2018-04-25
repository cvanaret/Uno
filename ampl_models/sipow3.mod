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

#   Source: problem 3 in
#   M. J. D. Powell,
#   "Log barrier methods for semi-infinite programming calculations"
#   Numerical Analysis Report DAMTP 1992/NA11, U. of Cambridge, UK.

#   SIF input: A. R. Conn and Nick Gould, August 1993

#   classification LLR2-AN-4-V

param M := 10000;
param STEP := 8/M;
param xi{j in 1..M/2} := if (1 <= j <= M/8) then 0.0
		    else if (M/8+1 <= j <= M/4) then (j-1)*STEP-1
		    else if (M/4+1 <= j <= 3*M/8) then 1
		    else if (3*M/8+1 <= j <= M/2) then (j-1)*STEP-3
		    ;
param eta{j in 1..M/2} := if (1 <= j <= M/8) then (j-1)*STEP
                    else if (M/8+1 <= j <= M/4) then 1
                    else if (M/4+1 <= j <= 3*M/8) then (j-1)*STEP-2
                    else if (3*M/8+1 <= j <= M/2) then 0
		    ;

param x_init{1..4};
var x{i in 1..4} := x_init[i];

minimize f:
	x[4];

subject to cons1{j in 1..M/2}:
	x[1]+x[4]+xi[j]*x[2]+eta[j]*x[3]-xi[j]^2*eta[j] >= 0;
subject to cons2{j in 1..M/2}:
	x[1]+xi[j]*x[2]+eta[j]*x[3]-xi[j]^2*eta[j] <= 0;

data;
param x_init:=
1	-0.1
2	0
3	0
4	1.2;

solve;
display f;
display x;
