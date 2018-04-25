param c := (48.4/50.176)*sin(.25);
param d := (48.4/50.176)*cos(.25);

var x {1..9};
var y1 = sin(x[8]);
var y2 = cos(x[8]);
var y3 = sin(x[9]);
var y4 = cos(x[9]);
var y5 = sin(x[8]-x[9]);
var y6 = cos(x[8]-x[9]);

minimize obj:
   3000*x[1]+1000*x[1]^3+2000*x[2]+666.667*x[2]^3;

s.t. c1:  .4-x[1]+2*c*x[5]^2-x[5]*x[6]*(d*y1+c*y2)-x[5]*x[7]*(d*y3+c*y4)=0;
s.t. c2:  .4-x[2]+2*c*x[6]^2+x[5]*x[6]*(d*y1-c*y2)+x[6]*x[7]*(d*y5-c*y6)=0;
s.t. c3:  .8+2*c*x[7]^2+x[5]*x[7]*(d*y3-c*y4)-x[6]*x[7]*(d*y5+c*y6)=0;
s.t. c4:  .2-x[3]+2*d*x[5]^2+x[5]*x[6]*(c*y1-d*y2)+x[5]*x[7]*(c*y3-d*y4)=0;
s.t. c5:  .2-x[4]+2*d*x[6]^2-x[5]*x[6]*(c*y1+d*y2)-x[6]*x[7]*(c*y5+d*y6)=0;
s.t. c6:  -.337+2*d*x[7]^2-x[5]*x[7]*(c*y3+d*y4)+x[6]*x[7]*(c*y5-d*y6)=0;
s.t. c7:  x[1]>=0;
s.t. c8:  x[2]>=0;
s.t. c9:  x[5]>=.90909;
s.t. c10: x[6]>=.90909;
s.t. c11: x[7]>=.90909;
s.t. c12: x[5]<=1.0909;
s.t. c13: x[6]<=1.0909;
s.t. c14: x[7]<=1.0909;

data;

let x[1] :=  .8;
let x[2] :=  .8;
let x[3] :=  .2;
let x[4] :=  .2;
let x[5] := 1.0454;
let x[6] := 1.0454;
let x[7] := 0;
let x[8] := 0;

display obj;

solve;

display x;

display obj;

display obj-5055.011803;



