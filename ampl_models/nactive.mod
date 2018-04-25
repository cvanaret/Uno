var x{1..2};

minimize f: x[1];
subject to c1: 0.5*(-x[1] - x[2]^2 - 1) >= 0;
subject to c2: x[1] - x[2]^2 >= 0;
subject to c3: -x[1] + x[2]^2 >= 0;

let x[1] := -20;
let x[2] := 10;
