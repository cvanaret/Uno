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
#   Ph. Toint, private communication.

#   SIF input: Ph. Toint, April 1991.

#   classification SOR2-AN-V-V

param t:=1000;
param xt := 10.0;
param mass := 0.37;
param tol := 0.1;

param h := xt/t;
param w := xt*(t+1)/2;

param fmax := xt/t;

var x{i in 0..t} >= 0.0, <= xt, := i*h;
var y{0..t};
var z{0..t};
var vx{0..t} := 1.0;
var vy{0..t};
var vz{0..t};

var ux{1..t} >= -fmax, <= fmax;
var uy{1..t} >= -fmax, <= fmax;
var uz{1..t} >= -fmax, <= fmax;

minimize f:
	sum {i in 1..t} (i*h/w)*(x[i] - xt)^2;
subject to acx{i in 1..t}:
	mass*(vx[i]-vx[i-1])/h - ux[i] = 0;
subject to acy{i in 1..t}:
        mass*(vy[i]-vy[i-1])/h - uy[i] = 0;
subject to acz{i in 1..t}:
        mass*(vz[i]-vz[i-1])/h - uz[i] = 0;
subject to psx{i in 1..t}:
        (x[i]-x[i-1])/h - vx[i] = 0;
subject to psy{i in 1..t}:
        (y[i]-y[i-1])/h - vy[i] = 0;
subject to psz{i in 1..t}:
        (z[i]-z[i-1])/h - vz[i] = 0;
subject to sc{i in 1..t}:
	(y[i] - sin(x[i]))^2 + (z[i] - cos(x[i]))^2 - tol^2 <= 0;

data;
fix x[0] := 0.0;
fix y[0] := 0.0;
fix z[0] := 1.0;
fix vx[0] := 0.0;
fix vy[0] := 0.0;
fix vz[0] := 0.0;
fix vx[t] := 0.0;
fix vy[t] := 0.0;
fix vz[t] := 0.0;

solve;
