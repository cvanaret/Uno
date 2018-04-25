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
#   F.S. Sisser,
#   "Elimination of bounds in optimization problems by transforming
#   variables",
#   Mathematical Programming 20:110-121, 1981.

#   See also Buckley#216 (p. 91)

#   SIF input: Ph. Toint, Dec 1989.

#   classification OUR2-AN-2-0

param xinit{1..2};
var x{i in 1..2} := xinit[i];

minimize f:
	3*x[1]^4 - 2*(x[1]*x[2])^2 + 3*x[2]^4;

data;
param xinit := 1 1.0 2 0.1;

solve; display f; display x;
