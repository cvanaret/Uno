#
#**************************
# SET UP THE INITIAL DATA *
#**************************
#   Problem :
#   *********
#   The problem is to find the minimum disc radius subject to polygon
#   determined by boundary discs intersecting all interior discs.
#   Source:
#   W. Pulleyblank,
#   private communication, 1991.
#   SIF input: A.R. Conn, November 1991.
#   classification LQR2-MY-29-23
#   Number of nodes
#   Number of sides to the polygon
#   Constants
#   Lower bound on the objective function
#   Objective function
#   Boundary constraints
#   Collinearity constraints
#ZL DISC2                              RNODES
#   Solution
param nnodes := 11;
param nlines := 6;
param x1 := 0.0;
param x2 := 8.0;
param x3 := 12.0;
param x4 := 8.0;
param x5 := 0.0;
param x6 := 4.0;
param x7 := 8.0;
param x8 := 8.0;
param x9 := 4.0;
param x10 := 2.0;
param x11 := 2.0;
param y1 := 10.0;
param y2 := 10.0;
param y3 := 5.0;
param y4 := 0.0;
param y5 := 0.0;
param y6 := 8.0;
param y7 := 7.0;
param y8 := 3.0;
param y9 := 1.0;
param y10 := 3.0;
param y11 := 6.0;
param rnodes := 11.0;

var epsilon >= 0.0 ,  <= 3.0, := 0.0;
var u1, := 0.0;
var v1, := 0.0;
var u2, := 0.0;
var v2, := 0.0;
var u3, := 0.0;
var v3, := 0.0;
var u4, := 0.0;
var v4, := 0.0;
var u5, := 0.0;
var v5, := 0.0;
var u6, := 0.0;
var v6, := 0.0;
var u7, := 0.0;
var v7, := 0.0;
var u8, := 0.0;
var v8, := 0.0;
var u9, := 0.0;
var v9, := 0.0;
var u10, := 0.0;
var v10, := 0.0;
var u11, := 0.0;
var v11, := 0.0;
var alpha1 >= 0.0 ,  <= 1.0, := 0.0;
var alpha2 >= 0.0 ,  <= 1.0, := 0.0;
var alpha3 >= 0.0 ,  <= 1.0, := 0.0;
var alpha4 >= 0.0 ,  <= 1.0, := 0.0;
var alpha5 >= 0.0 ,  <= 1.0, := 0.0;
var alpha6 >= 0.0 ,  <= 1.0, := 0.0;

minimize obj:
	epsilon;

subject to b1:
	(u1-0.0)^2 + (v1-10.0)^2 - epsilon^2 = 0;
subject to b2:
	(u2-8.0)^2 + (v2-10.0)^2 - epsilon^2 = 0;
subject to b3:
	(u3-12.0)^2 + (v3-5.0)^2 - epsilon^2 = 0;
subject to b4:
	(u4-8.0)^2 + (v4-0.0)^2 - epsilon^2 = 0;
subject to b5:
	(u5-0.0)^2 + (v5-0.0)^2 - epsilon^2 = 0;
subject to b6:
	0 >= (u6-4.0)^2 + (v6-8.0)^2 - epsilon^2;
subject to b7:
	0 >= (u7-8.0)^2 + (v7-7.0)^2 - epsilon^2;
subject to b8:
	0 >= (u8-8.0)^2 + (v8-3.0)^2 - epsilon^2;
subject to b9:
	0 >= (u9-4.0)^2 + (v9-1.0)^2 - epsilon^2;
subject to b10:
	0 >= (u10-2.0)^2 + (v10-3.0)^2 - epsilon^2;
subject to b11:
	0 >= (u11-2.0)^2 + (v11-6.0)^2 - epsilon^2;
subject to b162:
	(u2-u1) * (-alpha1) + u6 - u1 = 0;
subject to c162:
	(v2-v1) * (-alpha1) + v6 - v1 = 0;
subject to b273:
	(u3-u2) * (-alpha2) + u7 - u2 = 0;
subject to c273:
	(v3-v2) * (-alpha2) + v7 - v2 = 0;
subject to b384:
	(u4-u3) * (-alpha3) + u8 - u3 = 0;
subject to c384:
	(v4-v3) * (-alpha3) + v8 - v3 = 0;
subject to b495:
	(u5-u4) * (-alpha4) + u9 - u4 = 0;
subject to c495:
	(v5-v4) * (-alpha4) + v9 - v4 = 0;
subject to b5101:
	(u1-u5) * (-alpha5) + u10 - u5 = 0;
subject to c5101:
	(v1-v5) * (-alpha5) + v10 - v5 = 0;
subject to b5111:
	(u1-u5) * (-alpha5) + u11 - u5 = 0;
subject to c5111:
	(v1-v5) * (-alpha5) + v11 - v5 = 0;

solve;
display epsilon;
display u1;
display v1;
display u2;
display v2;
display u3;
display v3;
display u4;
display v4;
display u5;
display v5;
display u6;
display v6;
display u7;
display v7;
display u8;
display v8;
display u9;
display v9;
display u10;
display v10;
display u11;
display v11;
display alpha1;
display alpha2;
display alpha3;
display alpha4;
display alpha5;
display alpha6;
display obj;
