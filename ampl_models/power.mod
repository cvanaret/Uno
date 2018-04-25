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
#   S.S. Oren,
#   Self-scaling variable metric algorithms,
#   Part II: implementation and experiments"
#   Management Science 20(5):863-874, 1974.

#   See also Buckley#179 (p. 83)

#   SIF input: Ph. Toint, Dec 1989.

#   classification OUR2-AN-V-0

param N:=1000;
var x{1..N} := 1.0;

minimize f:
	sum {i in 1..N} (i*x[i])^2; 

solve; display f; display x;
