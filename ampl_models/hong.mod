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

#   Source: Se June Hong/Chid Apte

#   SIF input: A.R.Conn, Jan 1991.

#   classification OLR2-AN-4-1

var T{1..4} >= 0.0, <= 1.0, := 0.5;

minimize f:
	0.92+0.08*exp(0.38*25*T[1]) - 2.95+3.95*exp(0.11*50*T[2]) - 1.66+1657834*exp(-1.48*(9.0+4.0*T[3])) + 0.11+0.89*exp(0.00035*20000*T[4]);

subject to cons1:
	T[1]+T[2]+T[3]+T[4]-1.0 = 0;

solve; display f; display T;
