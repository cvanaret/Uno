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

#   Source: problem 2 in
#   M. J. D. Powell,
#   "Log barrier methods for semi-infinite programming calculations"
#   Numerical Analysis Report DAMTP 1992/NA11, U. of Cambridge, UK.

#   SIF input: A. R. Conn and Nick Gould, August 1993

#   classification LLR2-AN-2-V

param M:=10000;
param pi:=3.1415;
param xinit{1..2};
var x{i in 1..2}:=xinit[i];

minimize f:
	x[2];
subject to cons1{i in 1..M/2}:
	1+x[1]*cos(4*pi*i/M)+x[2]*sin(4*pi*i/M) >= 0;
subject to cons2:
	x[1]+1 >= 0;

data;
param xinit:= 1 0.8 2 0.5;

solve; display f; display x;
