param rho_ >= 0, default 1000;
param l{1..2} default 0;
param sl{1..2} default 0;
let l[1] := 700;
let sl[1] := 1;
var x {1..2};
var s{i in 1..2} >= sl[i];

minimize obj: 
	100*(x[2] - x[1]^2)^2 + (1 - x[1])^2
      - l[1]*(x[1]*x[2] - s[1])
      - l[2]*(x[1] + x[2]^2 - s[2])
      + rho_/2*( (x[1]*x[2] - s[1])^2 + (x[1] + x[2]^2 - s[2])^2 )
;

#subject to constr1: x[1]*x[2] - s[1] = 0;
#subject to constr2: x[1] + x[2]^2 - s[2] = 0;
subject to constr3: x[1] <= 1/2;

let x[1] := 0.45;
let x[2] :=  1.9;

solve;

display x,s,l;
display (x[1]*x[2] - s[1]);
display (x[1] + x[2]^2 - s[2]);

display obj;

display obj - 306.5;
