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

#   Source:  problem 8 in
#   A.R.Conn,N.I.M.Gould and Ph.L.Toint,
#   "Testing a class of methods for solving minimization 
#   problems with simple bounds on their variables, 
#   Mathematics of Computation 50, pp 399-430, 1988.

#   SIF input: Nick Gould and Ph. Toint, Dec 1995.

#   classification SUR2-AN-V-0

param ns := 499;
param n := 2*ns + 2;

var x{i in 1..n} := if (i>4) then -2.0;

minimize f:
	1.0 +
	sum {i in 1..ns} (
	100*(x[2*i]-x[2*i-1]^2)^2 +
	(1.0-x[2*i-1])^2 +
	90*(x[2*i+2]-x[2*i+1]^2)^2 +
	(1.0-x[2*i+1])^2 +
	10*(x[2*i]+x[2*i+2]-2.0)^2 +
	(x[2*i]-x[2*i+2])^2/10
	);

let x[1] := -3.0;
let x[2] := -1.0;
let x[3] := -3.0;
let x[4] := -1.0;

solve;
display f;
display x;
