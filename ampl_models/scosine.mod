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
#   N. Gould, private communication.

#   SIF input: N. Gould, Nov 1997

#   classification OUR2-AN-V-0

param N:=10000;
param scal := 12.0;
param scale{i in 1..N} := exp((i-1)*scal/(N-1));
var x{i in 1..N} := 1.0/scale[i];

minimize f:
	sum {i in 1..N-1}
cos(-0.5*scale[i+1]*x[i+1]+scale[i]^2*x[i]^2);

solve; display f; display x;
