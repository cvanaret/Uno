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

#   classification QBR2-AN-V-V

param N:=12;
param M:=6;
var x{1..N} >= 0.0, <= 10.0;
minimize f:
	sum {i in 1..N} -i*10*x[i] + sum {i in 1..M} x[i]*x[i+1];

solve;
