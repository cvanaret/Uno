var x {1..10} >= -100, <= 100;

param c {1..10};

minimize obj: 
  sum {j in 1..10} exp(x[j])*(c[j] + x[j] - log(sum {k in 1..10} exp(x[k])));

subject to cons1: 
  exp(x[1]) + 2*exp(x[2]) + 2*exp(x[3]) + exp(x[6]) + exp(x[10]) = 2;

subject to cons2: 
  exp(x[4]) + 2*exp(x[5]) + exp(x[6]) + exp(x[7]) = 1;

subject to cons3: 
  exp(x[3]) + exp(x[7]) + exp(x[8]) + 2*exp(x[9]) + exp(x[10]) = 1;

data;

param c := 
  1  -6.089
  2 -17.164
  3 -34.054
  4  -5.914
  5 -24.721
  6 -14.986
  7 -24.100
  8 -10.708
  9 -26.662
 10 -22.179
 ;

let {j in 1..10} x[j] := -2.3;

solve;

display x;

display obj;

display obj + 47.76109026;
