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
#   P.S. Ritch,
#   Automatica, 1973, V9, pp 415-429,
#   (example 6.1)

#   SIF input: Todd Plantenga, August 1995.

#   classification QQR2-AN-V-V

param t := 400;
param dt := 20/t;
param a1 := dt/10;
param springkm := 0.02;
param damping := 0.05;
param c1 := springkm*dt;
param c2 := damping*dt;

var x{i in 0..t} := 0.0;
var y{i in 0..t} >= -1.0, := -1.0;

var u{i in 0..t-1} >= -0.2, <= 0.2, := 0.0;

minimize f:
	a1*sum {i in 0..t} (x[i]^2) + 1000*y[t]^2;
subject to cons1{i in 0..t-1}:
	x[i+1]-x[i]-dt*y[i] = 0.0;
subject to cons2{i in 0..t-1}:
	y[i+1]-y[i]-dt*u[i]+c1*x[i]+c2*y[i]^3 = 0.0;

fix x[0] := 10.0;
fix y[0] := 0.0;
fix y[t] := 0.0;
solve;
display f;
display x,y,u;
