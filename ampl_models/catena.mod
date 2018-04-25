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
#   K. Veselic,
#   "De forma catenarum in campo gravitatis pendentium",
#   Klasicna Gimnazija u Zagrebu, Zagreb, 1987.

#   SIF input: Ph. L. Toint, May 1993.
#              correction by F. Ruediger, Mar 1997.

#   classification LQR2-AY-V-V

param N := 10;

param gamma := 9.81;
param tmass := 500.0;
param bl := 1.0;
param fract := 0.6;

param length := bl*(N+1)*fract;
param mass := tmass/(N+1);
param mg := mass*gamma;

var x{i in 0..N+1} := i*length/(N+1);
var y{i in 0..N+1} := -i*length/(N+1);
var z{0..N+1} := 0.0;

minimize f:
	mg*y[0]/2 + sum {i in 1..N} mg*y[i] + mg*y[N+1]/2;

subject to cons1{i in 1..N+1}:
	(x[i]-x[i-1])^2 + (y[i]-y[i-1])^2 + (z[i]-z[i-1])^2 = bl^2;

fix x[0] := 0.0;
fix y[0] := 0.0;
fix z[0] := 0.0;
fix x[N+1] := length;

solve;
display f;
display x,y,z;
