var x1 >= 1.0e-5, <= 2.0e+3;
var x2 >= 1.0e-5, <= 1.6e+4;
var x3 >= 1.0e-5, <= 1.2e+2;

var y {1..8};

param a{1..14};

minimize obj: -(0.063*y[2]*y[5] - 5.04*x1 - 3.36*y[3] - 0.035*x2 - 10*x3);

subject to constr1 {i in 1..7}:  y[i+1] >= a[i];
subject to constr2 {i in 8..14}: y[i-6] <= a[i];
subject to constr3: y[3] = 1.22*y[2] - x1;
subject to constr4: y[6] = (x2+y[3])/x1;
subject to constr5: y[2] = 0.01*x1*(112 + 13.167*y[6] - 0.6667*y[6]^2);
subject to constr6: y[5] = 86.35 + 1.098*y[6] - 0.038*y[6]^2 + 0.325*(y[4]-89);
subject to constr7: y[8] = 3*y[5] - 133;
subject to constr8: y[7] = 35.82 - 0.222*y[8];
subject to constr9: y[4] = 98000*x3/(y[2]*y[7] + 1000*x3);

data;

param a := 
    1	   0
    2	   0
    3	  85
    4	  90
    5	   3
    6	   0.01
    7	 145
    8	5000
    9	2000
    10	  93
    11	  95
    12	  12
    13	   4
    14	 162
    ;

let x1 := 1745;
let x2 := 12000;
let x3 := 110;

#let x1 := 1728.371286;
#let x2 := 16000;
#let x3 := 98.14151402;

param y2old;
let y[2] := 1.6*x1;
repeat {
    let y2old := y[2];
    let y[3] := 1.22*y[2] - x1;
    let y[6] := (x2+y[3])/x1;
    let y[2] := 0.01*x1*(112 + 13.167*y[6] - 0.6667*y[6]^2);
} while abs(y2old - y[2]) > 0.001 ;

param y4old;
let y[4] := 93;
repeat {
    let y4old := y[4];
    let y[5] := 86.35 + 1.098*y[6] - 0.038*y[6]^2 + 0.325*(y[4]-89);
    let y[8] := 3*y[5] - 133;
    let y[7] := 35.82 - 0.222*y[8];
    let y[4] := 98000*x3/(y[2]*y[7] + 1000*x3);
} while abs(y4old - y[4]) > 0.001 ;

display obj;

solve;

display obj;

display obj + 1162.036507;
