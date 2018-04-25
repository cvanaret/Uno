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
#   "The OPTIMA user manual (issue No.8, p. 37)",
#   Numerical Optimization Centre, Hatfield Polytechnic (UK), 1989.

#   SIF input: Ph. Toint, May 1990.

#   classification SUR2-AN-3-0

param t {1..21};
param z {1..21};

var x {1..3};

minimize obj: sum {j in 1..10} (exp(t[j]*x[3]) - x[1]*exp(t[j]*x[2]) + z[j])^2;

let x[1] :=  1.0;
let x[2] := -1.0;
let x[3] :=  0.0;

data;

param t :=
1 0.3
2 0.35
3 0.4
4 0.45
5 0.5
6 0.55
7 0.6
8 0.65
9 0.7
10 0.75
11 0.8
12 0.85
13 0.9
14 0.95
15 1.0
16 1.05
17 1.1
18 1.15
19 1.2
20 1.25
21 1.3
;

param z :=
1 1.561
2 1.473
3 1.391
4 1.313
5 1.239
6 1.169
7 1.103
8 1.04
9 0.981
10 0.925
11 0.8721
12 0.8221
13 0.7748
14 0.73
15 0.6877
16 0.6477
17 0.6099
18 0.5741
19 0.5403
20 0.5084
21 0.4782
;

solve;
display x;
display obj;
