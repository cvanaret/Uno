#
#**************************
# SET UP THE INITIAL DATA *
#**************************
#   Problem :
#   *********
#   A nonlinear optimal control problem with both state- and
#   control constraints.
#   The problem is to control (using an applied force of limited
#   magnitude) a mass moving in the plane, such that its
#   trajectory lies within a prescribed distance TOL of the curve y=sin(x),
#   and such that it stops at a given abscissa XT in minimum time.
#   The mass is initially stationary at the origin.
#   Source:
#   Ph. Toint, private communication.
#   SIF input: Ph. Toint, April 1991.
#   classification SOR2-AN-V-V
#   Number of time intervals
#   The number of variables is 6T+4, of which 6 are fixed.
#IE T                   50             $ n = 304
#IE T                   100            $ n = 604
#IE T                   500            $ n = 3004
#   Target abscissa
#   Mass
#   Tolerance along the sine trajectory
#   Constants
#   Useful parameters
#   Maximal force at any time
#   Solution
#LO SOLTN(10)           1.79999921
#LO SOLTN(50)           3.37869476
#LO SOLTN(100)          5.23053464
param t := 10;
param xt := 10.0;
param mass := 0.5;
param tol := 0.1;
param rt := 10.0;
param tp1 := 1.0 + (10.0);
param h := (10.0) / (10.0);
param mdh := (0.5) / ((10.0) / (10.0));
param xttp1 := (10.0) * (1.0 + (10.0));
param w := 0.5 * ((10.0) * (1.0 + (10.0)));
param ri := 10;
param ti := (10) * ((10.0) / (10.0));
param wdt1 := (0.5 * ((10.0) * (1.0 + (10.0)))) / ((1) * ((10.0) / (10.0)));
param wdt2 := (0.5 * ((10.0) * (1.0 + (10.0)))) / ((2) * ((10.0) / (10.0)));
param wdt3 := (0.5 * ((10.0) * (1.0 + (10.0)))) / ((3) * ((10.0) / (10.0)));
param wdt4 := (0.5 * ((10.0) * (1.0 + (10.0)))) / ((4) * ((10.0) / (10.0)));
param wdt5 := (0.5 * ((10.0) * (1.0 + (10.0)))) / ((5) * ((10.0) / (10.0)));
param wdt6 := (0.5 * ((10.0) * (1.0 + (10.0)))) / ((6) * ((10.0) / (10.0)));
param wdt7 := (0.5 * ((10.0) * (1.0 + (10.0)))) / ((7) * ((10.0) / (10.0)));
param wdt8 := (0.5 * ((10.0) * (1.0 + (10.0)))) / ((8) * ((10.0) / (10.0)));
param wdt9 := (0.5 * ((10.0) * (1.0 + (10.0)))) / ((9) * ((10.0) / (10.0)));
param wdt10 := (0.5 * ((10.0) * (1.0 + (10.0)))) / ((10) * ((10.0) / (10.0)));
param fmax := (10.0) / (10.0);
param im1 := -1 + (10);

var x0 >= 0.0 ,  <= 0.0;
var y0 >= 0.0 ,  <= 0.0;
var vx0 >= 0.0 ,  <= 0.0;
var vy0 >= 0.0 ,  <= 0.0;
var x1 >= 0.0 ,  <= 10.0 ,  := 1.0;
var y1;
var vx1 := 1.0;
var vy1;
var x2 >= 0.0 ,  <= 10.0 ,  := 2.0;
var y2;
var vx2 := 1.0;
var vy2;
var x3 >= 0.0 ,  <= 10.0 ,  := 3.0;
var y3;
var vx3 := 1.0;
var vy3;
var x4 >= 0.0 ,  <= 10.0 ,  := 4.0;
var y4;
var vx4 := 1.0;
var vy4;
var x5 >= 0.0 ,  <= 10.0 ,  := 5.0;
var y5;
var vx5 := 1.0;
var vy5;
var x6 >= 0.0 ,  <= 10.0 ,  := 6.0;
var y6;
var vx6 := 1.0;
var vy6;
var x7 >= 0.0 ,  <= 10.0 ,  := 7.0;
var y7;
var vx7 := 1.0;
var vy7;
var x8 >= 0.0 ,  <= 10.0 ,  := 8.0;
var y8;
var vx8 := 1.0;
var vy8;
var x9 >= 0.0 ,  <= 10.0 ,  := 9.0;
var y9;
var vx9 := 1.0;
var vy9;
var x10 >= 0.0 ,  <= 10.0 ,  := 10.0;
var y10;
var vx10 >= 0.0 ,  <= 0.0 ,  := 1.0;
var vy10 >= 0.0 ,  <= 0.0;
var ux1 >= -1.0 ,  <= 1.0;
var uy1 >= -1.0 ,  <= 1.0;
var ux2 >= -1.0 ,  <= 1.0;
var uy2 >= -1.0 ,  <= 1.0;
var ux3 >= -1.0 ,  <= 1.0;
var uy3 >= -1.0 ,  <= 1.0;
var ux4 >= -1.0 ,  <= 1.0;
var uy4 >= -1.0 ,  <= 1.0;
var ux5 >= -1.0 ,  <= 1.0;
var uy5 >= -1.0 ,  <= 1.0;
var ux6 >= -1.0 ,  <= 1.0;
var uy6 >= -1.0 ,  <= 1.0;
var ux7 >= -1.0 ,  <= 1.0;
var uy7 >= -1.0 ,  <= 1.0;
var ux8 >= -1.0 ,  <= 1.0;
var uy8 >= -1.0 ,  <= 1.0;
var ux9 >= -1.0 ,  <= 1.0;
var uy9 >= -1.0 ,  <= 1.0;
var ux10 >= -1.0 ,  <= 1.0;
var uy10 >= -1.0 ,  <= 1.0;

minimize obj:
	(((x1 - 10.0)*(x1 - 10.0))/55.0) + (((x2 - 10.0)*(x2 - 10.0))/27.5) + (((x3 - 10.0)*(x3 - 10.0))/18.333333333333332) + (((x4 - 10.0)*(x4 - 10.0))/13.75) + (((x5 - 10.0)*(x5 - 10.0))/11.0) + (((x6 - 10.0)*(x6 - 10.0))/9.166666666666666) + (((x7 - 10.0)*(x7 - 10.0))/7.857142857142857) + (((x8 - 10.0)*(x8 - 10.0))/6.875) + (((x9 - 10.0)*(x9 - 10.0))/6.111111111111111) + (((x10 - 10.0)*(x10 - 10.0))/5.5);

subject to acx1:
	0.5*vx1 - 0.5*vx0 - ux1 = 0;
subject to acy1:
	0.5*vy1 - 0.5*vy0 - uy1 = 0;
subject to psx1:
	x1 - x0 - vx1 = 0;
subject to psy1:
	y1 - y0 - vy1 = 0;
subject to sc1:
	0 <= -sin ( x1 )  + y1 + 0.1 <= 0.2;
subject to acx2:
	0.5*vx2 - 0.5*vx1 - ux2 = 0;
subject to acy2:
	0.5*vy2 - 0.5*vy1 - uy2 = 0;
subject to psx2:
	x2 - x1 - vx2 = 0;
subject to psy2:
	y2 - y1 - vy2 = 0;
subject to sc2:
	0 <= -sin ( x2 )  + y2 + 0.1 <= 0.2;
subject to acx3:
	0.5*vx3 - 0.5*vx2 - ux3 = 0;
subject to acy3:
	0.5*vy3 - 0.5*vy2 - uy3 = 0;
subject to psx3:
	x3 - x2 - vx3 = 0;
subject to psy3:
	y3 - y2 - vy3 = 0;
subject to sc3:
	0 <= -sin ( x3 )  + y3 + 0.1 <= 0.2;
subject to acx4:
	0.5*vx4 - 0.5*vx3 - ux4 = 0;
subject to acy4:
	0.5*vy4 - 0.5*vy3 - uy4 = 0;
subject to psx4:
	x4 - x3 - vx4 = 0;
subject to psy4:
	y4 - y3 - vy4 = 0;
subject to sc4:
	0 <= -sin ( x4 )  + y4 + 0.1 <= 0.2;
subject to acx5:
	0.5*vx5 - 0.5*vx4 - ux5 = 0;
subject to acy5:
	0.5*vy5 - 0.5*vy4 - uy5 = 0;
subject to psx5:
	x5 - x4 - vx5 = 0;
subject to psy5:
	y5 - y4 - vy5 = 0;
subject to sc5:
	0 <= -sin ( x5 )  + y5 + 0.1 <= 0.2;
subject to acx6:
	0.5*vx6 - 0.5*vx5 - ux6 = 0;
subject to acy6:
	0.5*vy6 - 0.5*vy5 - uy6 = 0;
subject to psx6:
	x6 - x5 - vx6 = 0;
subject to psy6:
	y6 - y5 - vy6 = 0;
subject to sc6:
	0 <= -sin ( x6 )  + y6 + 0.1 <= 0.2;
subject to acx7:
	0.5*vx7 - 0.5*vx6 - ux7 = 0;
subject to acy7:
	0.5*vy7 - 0.5*vy6 - uy7 = 0;
subject to psx7:
	x7 - x6 - vx7 = 0;
subject to psy7:
	y7 - y6 - vy7 = 0;
subject to sc7:
	0 <= -sin ( x7 )  + y7 + 0.1 <= 0.2;
subject to acx8:
	0.5*vx8 - 0.5*vx7 - ux8 = 0;
subject to acy8:
	0.5*vy8 - 0.5*vy7 - uy8 = 0;
subject to psx8:
	x8 - x7 - vx8 = 0;
subject to psy8:
	y8 - y7 - vy8 = 0;
subject to sc8:
	0 <= -sin ( x8 )  + y8 + 0.1 <= 0.2;
subject to acx9:
	0.5*vx9 - 0.5*vx8 - ux9 = 0;
subject to acy9:
	0.5*vy9 - 0.5*vy8 - uy9 = 0;
subject to psx9:
	x9 - x8 - vx9 = 0;
subject to psy9:
	y9 - y8 - vy9 = 0;
subject to sc9:
	0 <= -sin ( x9 )  + y9 + 0.1 <= 0.2;
subject to acx10:
	0.5*vx10 - 0.5*vx9 - ux10 = 0;
subject to acy10:
	0.5*vy10 - 0.5*vy9 - uy10 = 0;
subject to psx10:
	x10 - x9 - vx10 = 0;
subject to psy10:
	y10 - y9 - vy10 = 0;
subject to sc10:
	0 <= -sin ( x10 )  + y10 + 0.1 <= 0.2;

solve;
display x0;
display y0;
display vx0;
display vy0;
display x1;
display y1;
display vx1;
display vy1;
display x2;
display y2;
display vx2;
display vy2;
display x3;
display y3;
display vx3;
display vy3;
display x4;
display y4;
display vx4;
display vy4;
display x5;
display y5;
display vx5;
display vy5;
display x6;
display y6;
display vx6;
display vy6;
display x7;
display y7;
display vx7;
display vy7;
display x8;
display y8;
display vx8;
display vy8;
display x9;
display y9;
display vx9;
display vy9;
display x10;
display y10;
display vx10;
display vy10;
display ux1;
display uy1;
display ux2;
display uy2;
display ux3;
display uy3;
display ux4;
display uy4;
display ux5;
display uy5;
display ux6;
display uy6;
display ux7;
display uy7;
display ux8;
display uy8;
display ux9;
display uy9;
display ux10;
display uy10;
display obj;
