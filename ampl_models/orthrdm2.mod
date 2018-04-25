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

#   Source: adapted from:
#   M. Gulliksson,
#   "Algorithms for nonlinear Least-squares with Applications to
#   Orthogonal Regression",
#   UMINF-178.90, University of Umea, Sweden, 1990.

#   SIF input: Ph. Toint, Mar 1991,
#              modified by T. Plantagena, May 1994.

#   classification QOR2-AY-V-V

param npts := 2000;
param tz3 := 1.7;
param pseed := 237.1531;
param psize := 0.1;
param pi := 4*atan(1);
param icr0 := 1/npts;
param incr := icr0*2*pi;

param xd{i in 1..npts} := ((1+tz3^2)+cos(incr*(i-1)))*cos(incr*(i-1))*(1+psize*cos(incr*(i-1)*pseed));
param yd{i in 1..npts} := ((1+tz3^2)+cos(incr*(i-1)))*sin(incr*(i-1))*(1+psize*cos(incr*(i-1)*pseed));

var z1 := 1.0;
var z2 := 0.0;
var z3 := 1.0;
var x{i in 1..npts} := xd[i];
var y{i in 1..npts} := yd[i];

minimize f:
	sum {i in 1..npts} (
	(x[i]-xd[i])^2 +
	(y[i]-yd[i])^2 );
subject to cons1{i in 1..npts}:
	((x[i]-z1)^2+(y[i]-z2)^2)^2 - 
	((x[i]-z1)^2+(y[i]-z2)^2)*(1+z3^2)^2= 0.0;

solve;
display f;
display x, y, z1, z2, z3;

