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

#   Source: An inverse problem from astronomy,
#   reformulated as a convex quadratic program by
#   S. P. Hestis, SIAM Review 34 (1992) pp. 642-647.

#   SIF input: Nick Gould, January 1993.
#   improvements by: Ruediger Franke (Ruediger.Franke@RZ.TU-Ilmenau.DE)

#   classification QLR2-MN-V-V

param K := 10000;
param range := 1.0;
param deltax := range/K;
var M{1..K} >= 0, :=1.0;

minimize f:
	sum {i in 1..K} M[i]^2/K;
subject to cons1:
	sum {i in 1..K} ((i^3)-((i-1)^3))*(deltax^3)*M[i]/3 - 1835.2 = 0;
subject to cons2:
	sum {i in 1..K} ((i^5)-((i-1)^5))*(deltax^5)*M[i]/5 - 909.8 = 0;

solve;
display f;
display M;
