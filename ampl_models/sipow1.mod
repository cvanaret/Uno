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

#   Source: problem 1 in
#   M. J. D. Powell,
#   "Log barrier methods for semi-infinite programming calculations"
#   Numerical Analysis Report DAMTP 1992/NA11, U. of Cambridge, UK.

#   SIF input: A. R. Conn and Nick Gould, August 1993

#   classification LLR2-AN-2-V

param m := 10000;
param xinit{1..2};
param pi := 3.1415;
var x{i in 1..2} := xinit[i];

minimize f:
	x[2];
subject to cons{j in 1..m}:
	x[1]*cos(2*pi*j/m) + x[2]*sin(2*pi*j/m) + 1.0 >= 0;

data;
param xinit:= 1 0.8 2 0.5;

solve;
display f;
display x;
