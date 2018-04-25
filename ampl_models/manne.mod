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
#   B. Murtagh and M. Saunders,
#   Mathematical Programming Studies 16, pp. 84-117,
#   (example 5.12).

#   SIF input: N. Gould and Ph. Toint, March 1990.

#   classification OOR2-MN-V-V

param t :=365;

param grow := 0.03;
param beta := 0.95;
param xk0 := 3.0;
param xc0 := 0.95;
param xi0 := 0.05;
param b := 0.25;
param bprob := 1.0;

param a := (xc0 + xi0)/(xk0^b);
param gfac := (1+grow)^(1-b);

param at{j in 1..t} := if (j == 1) then a*gfac else at[j-1]*gfac;
param bt{j in 1..t} := if (j == 1) then beta else if (2 <= j <= t-1) then  bt[j-1]*beta else bt[t-1]*beta/(1-beta);

var c{1..t} >= 0.95, := 0.95;
var ii{i in 1..t} >= 0.05, <= 0.05*(1.04)^i, := 0.05;
var kk{i in 1..t} >= 3.05, := if (2 <= i <= t) then (3+(i-1)/10);

minimize f:
	sum {i in 1..t} bt[i]*log(c[i]);
subject to cons1{i in 1..t}:
	at[i]*kk[i]^b- c[i] - ii[i] >= 0;
subject to cons2{i in 1..t-1}:
	kk[i+1] - kk[i] - ii[i] <= 0;
subject to cons3:
	grow*kk[t] - ii[t] <= 0;

fix kk[1] := 3.05;
solve;
	
