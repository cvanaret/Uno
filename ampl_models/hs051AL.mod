var x {1..5};
param l{1..3} default 0;
param rho_ := 100;

minimize obj: 
  (x[1]-x[2])^2 + (x[2]+x[3]-2)^2 + (x[4]-1)^2 + (x[5]-1)^2
  - l[1]*(x[1] + 3*x[2] - 4)
  - l[2]*(x[3] +   x[4] - 2*x[5])
  - l[3]*(x[2] -   x[5])
  + rho_*( (x[1] + 3*x[2] - 4)^2 + (x[3] +   x[4] - 2*x[5])^2 + (x[2] -   x[5])^2 )
;

let x[1] :=  2.5;
let x[2] :=  0.5;
let x[3] :=  2;
let x[4] := -1;
let x[5] :=  0.5;

display obj;

solve;

display x;

display obj, rho_*( (x[1] + 3*x[2] - 4)^2 + (x[3] +   x[4] - 2*x[5])^2 + (x[2] -   x[5])^2 );

display obj - 0;
