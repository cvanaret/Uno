var x{1..2};

minimize f: x[1] + x[2];
subject to c1: x[2] - x[1]^2 - 1 >= 0;
subject to c2: 0.3*(1 - exp(x[2])) >= 0;

let x[1] := 3;
let x[2] := 2;
