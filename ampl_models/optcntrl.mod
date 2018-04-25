#
#**************************
# SET UP THE INITIAL DATA *
#**************************
#   Problem :
#   *********
#   An optimal control problem
#   Source:
#   B. Murtagh and M. Saunders,
#   Mathematical Programming studies 16, pp 84-117,
#   (example 5.11)
#   SIF input: Nick Gould, June 1990.
#   classification QQR2-AN-32-20
#   useful parameters
#   Least square problems are bounded below by zero
#   Solution
param t := 10;
param tm1 := -1 + (10);
param tp1 := 1 + (9);

var x0 >= 10.0 ,  <= 10.0;
var y0 >= 0.0 ,  <= 0.0;
var x1;
var y1 >= -1.0;
var x2;
var y2 >= -1.0;
var x3;
var y3 >= -1.0;
var x4;
var y4 >= -1.0;
var x5;
var y5 >= -1.0;
var x6;
var y6 >= -1.0;
var x7;
var y7 >= -1.0;
var x8;
var y8 >= -1.0;
var x9;
var y9 >= -1.0;
var x10 >= 0.0;
var y10 >= 0.0 ,  <= 0.0;
var u0 >= -0.2 ,  <= 0.2;
var u1 >= -0.2 ,  <= 0.2;
var u2 >= -0.2 ,  <= 0.2;
var u3 >= -0.2 ,  <= 0.2;
var u4 >= -0.2 ,  <= 0.2;
var u5 >= -0.2 ,  <= 0.2;
var u6 >= -0.2 ,  <= 0.2;
var u7 >= -0.2 ,  <= 0.2;
var u8 >= -0.2 ,  <= 0.2;
var u9 >= -0.2 ,  <= 0.2;

minimize obj:
	0.5*x0 * x0 + 0.5*x1 * x1 + 0.5*x2 * x2 + 0.5*x3 * x3 + 0.5*x4 * x4 + 0.5*x5 * x5 + 0.5*x6 * x6 + 0.5*x7 * x7 + 0.5*x8 * x8 + 0.5*x9 * x9 + 0.5*x10 * x10;

subject to obj_bnd:
	0.0 <= 0.5*x0 * x0 + 0.5*x1 * x1 + 0.5*x2 * x2 + 0.5*x3 * x3 + 0.5*x4 * x4 + 0.5*x5 * x5 + 0.5*x6 * x6 + 0.5*x7 * x7 + 0.5*x8 * x8 + 0.5*x9 * x9 + 0.5*x10 * x10;
subject to b0:
	x1 - x0 - 0.2*y0 = 0;
subject to c0:
	0.01*y0 * y0 + y1 - y0 + 0.0040*x0 - 0.2*u0 = 0;
subject to b1:
	x2 - x1 - 0.2*y1 = 0;
subject to c1:
	0.01*y1 * y1 + y2 - y1 + 0.0040*x1 - 0.2*u1 = 0;
subject to b2:
	x3 - x2 - 0.2*y2 = 0;
subject to c2:
	0.01*y2 * y2 + y3 - y2 + 0.0040*x2 - 0.2*u2 = 0;
subject to b3:
	x4 - x3 - 0.2*y3 = 0;
subject to c3:
	0.01*y3 * y3 + y4 - y3 + 0.0040*x3 - 0.2*u3 = 0;
subject to b4:
	x5 - x4 - 0.2*y4 = 0;
subject to c4:
	0.01*y4 * y4 + y5 - y4 + 0.0040*x4 - 0.2*u4 = 0;
subject to b5:
	x6 - x5 - 0.2*y5 = 0;
subject to c5:
	0.01*y5 * y5 + y6 - y5 + 0.0040*x5 - 0.2*u5 = 0;
subject to b6:
	x7 - x6 - 0.2*y6 = 0;
subject to c6:
	0.01*y6 * y6 + y7 - y6 + 0.0040*x6 - 0.2*u6 = 0;
subject to b7:
	x8 - x7 - 0.2*y7 = 0;
subject to c7:
	0.01*y7 * y7 + y8 - y7 + 0.0040*x7 - 0.2*u7 = 0;
subject to b8:
	x9 - x8 - 0.2*y8 = 0;
subject to c8:
	0.01*y8 * y8 + y9 - y8 + 0.0040*x8 - 0.2*u8 = 0;
subject to b9:
	x10 - x9 - 0.2*y9 = 0;
subject to c9:
	0.01*y9 * y9 + y10 - y9 + 0.0040*x9 - 0.2*u9 = 0;

solve;
display x0;
display y0;
display x1;
display y1;
display x2;
display y2;
display x3;
display y3;
display x4;
display y4;
display x5;
display y5;
display x6;
display y6;
display x7;
display y7;
display x8;
display y8;
display x9;
display y9;
display x10;
display y10;
display u0;
display u1;
display u2;
display u3;
display u4;
display u5;
display u6;
display u7;
display u8;
display u9;
display obj;
