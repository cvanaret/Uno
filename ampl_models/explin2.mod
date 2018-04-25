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

#   classification OBR2-AN-V-V

param n := 120;
param m := 10;

var x{1..n} >= 0, <= 10.0, := 0.0;

minimize f:
	sum {i in 1..m} exp(0.1*i*x[i]*x[i+1]/m) + sum {i in 1..n} (-10.0*i*x[i]);

solve;

