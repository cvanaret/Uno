###############################################################
# I can't get the initial objective function value(s) to match
# those reported in Hock-Schittkowski.  I don't see why.
###############################################################

param c {1..19};
param u {1..4};
param y_obs {1..19};

var x {j in 1..4} >= 0.00001, <= u[j];
var b = x[3] + (1-x[3])*x[4];
var y_cal {i in 1..19} = 
    (1 + 1/(12*x[2]))
    *
    (
	x[3]*b^x[2]*(x[2]/6.2832)^0.5 * (c[i]/7.685)^(x[2]-1)
	* exp(x[2] - b*c[i]*x[2]/7.658)
    )
    +
    (1 + 1/(12*x[1]))
    *
    (
	(1-x[3])*(b/x[4])^x[1]*(x[1]/6.2832)^0.5 * (c[i]/7.658)^(x[1]-1)
	* exp(x[1] - b*c[i]*x[1]/(7.658*x[4]))
    );

minimize obj: 
  sum {i in 1..19} (y_cal[i] - y_obs[i])^2
  ;

subject to constr1: x[3] + (1-x[3])*x[4] >= 0;

data;

param u :=
 1  100
 2  100
 3    1
 4  100
 ;

param: c y_obs :=
 1 .1 0.00189
 2  1 0.1038
 3  2 0.268
 4  3 0.506
 5  4 0.577
 6  5 0.604
 7  6 0.725
 8  7 0.898
 9  8 0.947
10  9 0.845
11 10 0.702
12 11 0.528
13 12 0.385
14 13 0.257
15 14 0.159
16 15 0.0869
17 16 0.0453
18 17 0.01509
19 18 0.00189
;

let x[1] := 2;
let x[2] := 4;
let x[3] := 0.04;
let x[4] := 2;

#printf "optimal solution as starting point \n";
#let x[1] := 12.27695;
#let x[2] := 4.631788;
#let x[3] := 0.3128625;
#let x[4] := 2.029290;

display obj;

solve;

display x;

display obj;

display obj - 0.007498464;
