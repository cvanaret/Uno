param rho_ >= 0, default 1000;
param l{1..3} default 0;
param sl{1..3} default 0;
let l[1] := 0.9;
let l[3] := 0.5;
let sl[3] := -6;

var x {1..2};
var s {i in 1..3} >= sl[i];


minimize obj: 
  ((x[1] - 3)^2 - 9) * x[2]^3 / (27*sqrt(3))
-l[1]*( x[1]/sqrt(3) - x[2] - s[1] )
-l[2]*( x[1] + sqrt(3)*x[2] - s[2] )
-l[3]*(-x[1] - sqrt(3)*x[2] - s[3] )
+rho_/2*( (x[1]/sqrt(3) - x[2] - s[1] )^2 + (x[1] + sqrt(3)*x[2] - s[2] )^2 + (-x[1] - sqrt(3)*x[2] - s[3])^2 )
;

let x[1] := 1;
let x[2] := 1/2;

#printf "optimal solution as starting point \n";
#let x[1] := 3;
#let x[2] := 1.73205;





let x[1] := 2;
let x[2] := 2;

#printf "optimal solution as starting point \n";
#let x[1] := 0.822876;
#let x[2] := 0.911438;

solve;

display x,s,l;

display obj;

display obj - 9 + 2.875*sqrt(7);
