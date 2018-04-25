param a {1..44};
param b {1..44};

var x {1..2};

minimize obj: 
  sum {i in 1..44} (b[i] - x[1] - (0.49 - x[1])*exp(-x[2]*(a[i]-8)))^2
  ;

subject to constr1: 0.49*x[2] - x[1]*x[2] >= 0.09;
subject to constr2: x[1] >= 0.4;
subject to constr3: x[2] >= -4;

let x[1] := 0.42;
let x[2] := 5;

#printf "optimal solution as starting point \n";
#let x[1] := 0.419952675;
#let x[2] := 1.284845629;

data;
param: a b :=
 1  8 0.49
 2  8 0.49
 3 10 0.48
 4 10 0.47
 5 10 0.48
 6 10 0.47
 7 12 0.46
 8 12 0.46
 9 12 0.45
10 12 0.43
11 14 0.45
12 14 0.43
13 14 0.43
14 16 0.44
15 16 0.43
16 16 0.43
17 18 0.46
18 18 0.45
19 20 0.42
20 20 0.42
21 20 0.43
22 22 0.41
23 22 0.41
24 22 0.40
25 24 0.42
26 24 0.40
27 24 0.40
28 26 0.41
29 26 0.40
30 26 0.41
31 28 0.41
32 28 0.40
33 30 0.40
34 30 0.40
35 30 0.38
36 32 0.41
37 32 0.40
38 34 0.40
39 36 0.41
40 36 0.38
41 38 0.40
42 38 0.40
43 40 0.39
44 42 0.39
;

/*
display obj;
display constr1.lb, constr1.body, constr1.ub;
display constr2.lb, constr2.body, constr2.ub;
display constr3.lb, constr3.body, constr3.ub;
*/

solve;

display x;

display obj;

display obj - 0.02845966972;
