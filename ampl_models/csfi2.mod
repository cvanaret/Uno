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

#   Source: problem MINLEN in
#   Vasko and Stott
#   "Optimizing continuous caster product dimensions:
#    an example of a nonlinear design problem in the steel industry"
#   SIAM Review, Vol 37 No, 1 pp.82-84, 1995

#   SIF input: A.R. Conn April 1995

#   classification LOR2-RN-5-4

param mintph := 45.0;
param minthick := 7.0;
param minarea := 200.0;
param maxarea := 250.0;
param maxaspr := 2.0;
param k := 1.0;

var thick >= minthick, := 0.5;
var wid >= 0.0, := 0.5;
var len >= 0.0, := 0.5;
var tph >= mintph, := 0.5;
var ipm >= 0.0, := 0.5;

minimize f:
	len;
subject to cons1:
	117.370892*tph/(wid*thick)-ipm = 0.0;
subject to cons2:
	thick^2*ipm/48.0-len = 0.0;
subject to cons3:
	wid/thick <= maxaspr;
subject to cons4:
	0.0 <= thick*wid - minarea <= maxarea-minarea;

solve;
display f;
display thick, wid, len, tph, ipm;
