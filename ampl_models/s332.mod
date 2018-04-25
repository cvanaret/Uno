param N := 2;
param pi := 4*atan(1); # 3.1415;
var x{1..N};
param t{i in 1..100} := pi*(1/3+(i-1)/180);
#var pmax = max {i in 1..100} ((180/pi)*atan((1/t[i]-x[1])/(log(t[i])+x[2])));

minimize f:
(pi/3.6)*sum {i in 1..100} ((log(t[i])+x[2]*sin(t[i])+x[1]*cos(t[i]))^2 +
(log(t[i])+x[2]*cos(t[i])-x[1]*sin(t[i]))^2);

subject to cons1 {i in 1..100}:
-30 <= (180/pi)*atan((1/t[i]-x[1])/(log(t[i])+x[2])) <= 30;

subject to cons2{i in 1..2}:
0 <= x[i] <= 1.5;

data;
var x:=
1	0.75
2	0.75;

solve;

display x;
