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

#   Source:  adapted from:
#   M. Gulliksson,
#   "Algorithms for nonlinear Least-squares with Applications to
#   Orthogonal Regression",
#   UMINF-178.90, University of Umea, Sweden, 1990.

#   SIF input: Ph. Toint, June 1990.

#   classification QQR2-AN-V-V

param npts := 5000;
param v11 := 2.0;
param v12 := 1.0;
param v22 := 2.0;

param pseed := 237.1531;
param psize := 0.2;

param pi := 4*atan(1);

param incr := 2*pi/npts;
param c3 := cos(v22);
param s3 := sin(v22);

param xd{i in 1..npts} := (v11*cos(incr*(i-1))*c3-v12*sin(incr*(i-1))*s3)*(psize*cos(incr*(i-1)*pseed)+1);
param yd{i in 1..npts} := (v11*cos(incr*(i-1))*s3+v12*sin(incr*(i-1))*c3)*(psize*cos(incr*(i-1)*pseed)+1);

var h11 := 1.0;
var h12;
var h22 := 1.0;
var g1 := 1.0;
var g2 := 1.0;

var x{i in 1..npts} := xd[i];
var y{i in 1..npts} := yd[i];

minimize f:
	sum {i in 1..npts} (x[i]-xd[i])^2 
	+ sum {i in 1..npts} (y[i]-yd[i])^2;

subject to cons1{i in 1..npts}:
	h11*x[i]^2 + 2*h12*x[i]*y[i] + h22*y[i]^2 - 2*g1*x[i] - 2*g2*y[i] = 1.0;

solve;
display f;
display h11,h12,h22,g1,g2;
display x,y;
