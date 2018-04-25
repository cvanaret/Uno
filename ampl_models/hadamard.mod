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

#   Source:  A suggestion by Alan Edelman (MIT).

#   SIF input: Nick Gould, Nov 1993.

#   classification LQR2-RN-V-V

param n := 8;

set N := 1..n;

var Q{N,N} := 1;

var maxval >= 0, := 0;

var QQT{i in N, j in N} = sum{k in N} (Q[k,i]*Q[k,j]);

minimize abs_val:
	abs(maxval);

subject to ortho{i in N, j in N}:
	QQT[i,j] - n = 0;

subject to abs_min_val{i in N, j in N}:
	maxval >= Q[i,j];

subject to abs_max_val{i in N, j in N}:
	maxval >= -Q[i,j];

subject to ones{i in N,j in N}:
	abs(Q[i,j]) <= 1;

solve;  

display Q;
