var x {1..10} >= 1.0e-6;

param c {1..10};

minimize obj: 
  sum {j in 1..10} x[j]*(c[j] 
    + log(x[j]/sum {k in 1..10} x[k]) );

subject to cons1: 
  x[1] + 2*x[2] + 2*x[3] + x[6] + x[10] = 2;

subject to cons2: 
  x[4] + 2*x[5] + x[6] + x[7] = 1;

subject to cons3: 
  x[3] + x[7] + x[8] + 2*x[9] + x[10] = 1;

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

let {j in 1..10} x[j] := 0.1;

solve;

display x;

display obj;

display obj + 47.76109026;
