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

#   Source: problem 100 in
#   W. Hock and K. Schittkowski,
#   "Test examples for nonlinear programming codes",
#   Lectures Notes in Economics and Mathematical Systems 187, Springer
#   Verlag, Heidelberg, 1981.
#   This problem has been modified 20 Oct 92 by Todd Plantenga as follows.
#   The nonlinear inequality constraints are removed (if inactive
#   at the solution) or changed to equalities (if active).

#   SIF input: Ph. Toint, April 1991 and T. Plantenga, October 1992.

#   classification OOR2-AN-7-2

var x {1..7};

minimize obj: 
  (x[1]-10)^2 + 5*(x[2]-12)^2 + x[3]^4 + 3*(x[4]-11)^2 + 10*x[5]^6
  + 7*x[6]^2 + x[7]^4 - 4*x[6]*x[7] - 10*x[6] - 8*x[7]
  ;

subject to constr1: 2*x[1]^2 + 3*x[2]^4 + x[3] + 4*x[4]^2 + 5*x[5] = 127;
subject to constr4: -4*x[1]^2 - x[2]^2 + 3*x[1]*x[2] -2*x[3]^2 - 5*x[6]
			+11*x[7] = 0;

data;

let x[1] := 1;
let x[2] := 2;
let x[3] := 0;
let x[4] := 4;
let x[5] := 0;
let x[6] := 1;
let x[7] := 1;

#printf "optimal solution as starting point \n";
#let x[1] := 2.330499;
#let x[2] := 1.951372;
#let x[3] := -0.4775414;
#let x[4] := 4.365726;
#let x[5] := -0.6244870;
#let x[6] := 1.038131;
#let x[7] := 1.594227;


solve;

display x;

display obj;

display obj - 680.6300573;
