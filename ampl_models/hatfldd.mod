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
#   "The OPTIMA user manual (issue No.8, p. 35)",
#   Numerical Optimization Centre, Hatfield Polytechnic (UK), 1989.

#   SIF input: Ph. Toint, May 1990.

#   classification SUR2-AN-3-0

param t {1..10};
param z {1..10};

var x {1..3};

minimize obj: sum {j in 1..10} (exp(t[j]*x[3]) - x[1]*exp(t[j]*x[2]) + z[j])^2;

let x[1] :=  1.0;
let x[2] := -1.0;
let x[3] :=  0.0;

data;
param t :=
 1 0.2
 2 0.3
 3 0.4
 4 0.5
 5 0.6
 6 0.7
 7 0.75
 8 0.8
 9 0.85
 10 0.9 ;

param z :=
1 1.751
2 1.561
3 1.391
4 1.239
5 1.103
6 0.981
7 0.925
8 0.8721
9 0.8221
10 0.7748 ;

solve;
display x;
display obj;
