var x {1..3} >=0, <= 1;


minimize obj: 
  -32.174*(255*log((x[1]+x[2]+x[3]+0.03)/(0.09*x[1] + x[2] + x[3] + 0.03))
	  +280*log((x[2]+x[3]+0.03)/(0.07*x[2] + x[3] + 0.03))
	  +290*log((x[3]+0.03)/(0.13*x[3] + 0.03)))
  ;

subject to constr1: x[1] + x[2] + x[3] = 1;

let x[1] := 0.7;
let x[2] := 0.2;
let x[3] := 0.1;

#printf "optimal solution as starting point \n";
#let x[1] := 0.6178126908;
#let x[2] := 0.328202223;
#let x[3] := 0.5398508606e-1;

display obj;

solve;

display x;

display obj;

display obj + 26272.51448;
