#
#**************************
# SET UP THE INITIAL DATA *
#**************************
#   Problem :
#   *********
#   An optimal control version of the Singly Supported NonLinear BEAM problem.
#   The energy of a beam of length 1 compressed by a force P is to be
#   minimized.  The control variable is the derivative of the deflection angle.
#   The problem is discretized using the trapezoidal rule. It is non convex.
#   Source:
#   H. Maurer and H.D. Mittelman,
#   "The non-linear beam via optimal control with bound state variables",
#   Optimal Control Applications and Methods 12, pp. 19-31, 1991.
#   SIF input: Ph. Toint, Nov 1993.
#   classification  OOR2-MN-V-V
#   Discretization: specify the number of interior points + 1
#IE NI                  50             $ n =   153, m =   100
#IE NI                  100            $ n =   303, m =   200
#IE NI                  500            $ n =  1503, m =  1000
#IE NI                  1000           $ n =  3003, m =  1000
#IE NI                  5000           $ n = 15003, m = 10000
#  Set ALPHA to the force divided by the bending stiffness.
#  For this value, 3 contact zones with the boundary of maximum displacement
#  occur, alternatively with negative and positive displacement.
#   Useful constants
#   deflection angle
#   displacement
#   derivative of the deflection angle (control variable)
#   objective function: energy
#   state equations
#   All variables are free
#   Impose the bound on the displacements
#   Impose that the beam does not curve too far backwards
#XL SSNLBEAM  T(I)      -1.57079663
#XU SSNLBEAM  T(I)      1.57079663
#   Boundary conditions
#  The origin is a saddle point!
#  Perturb the origin
#   Solution
#LO SOLTN(10)           337.7716271268
#LO SOLTN(50)           340.0150353685
#LO SOLTN(100)          340.0284918253
#LO SOLTN(500)          340.0323978636
#LO SOLTN(1000)         
	param ni := 10;
	param alpha := 350.0;
	param rni := 10.0;
	param nim1 := -1 + (10);
	param h := 1.0 / (10.0);
	param hd4 := 0.25 * (1.0 / (10.0));
	param hd2 := 0.5 * (1.0 / (10.0));
	param ah := (350.0) * (1.0 / (10.0));
	param ahd2 := 0.5 * ((350.0) * (1.0 / (10.0)));
	param ip1 := 1 + (9);
	param ri := 10;
	param tt := (10) * (1.0 / (10.0));
	param ctt := cos((10) * (1.0 / (10.0)));
	param sctt := 0.05 * (cos((10) * (1.0 / (10.0))));

	var t0 >= -1.0 ,  <= 1.0 ,  := 0.05;
	var t1 >= -1.0 ,  <= 1.0 ,  := 0.04975020826390129;
	var t2 >= -1.0 ,  <= 1.0 ,  := 0.04900332889206208;
	var t3 >= -1.0 ,  <= 1.0 ,  := 0.047766824456280305;
	var t4 >= -1.0 ,  <= 1.0 ,  := 0.04605304970014426;
	var t5 >= -1.0 ,  <= 1.0 ,  := 0.04387912809451864;
	var t6 >= -1.0 ,  <= 1.0 ,  := 0.041266780745483914;
	var t7 >= -1.0 ,  <= 1.0 ,  := 0.038242109364224425;
	var t8 >= -1.0 ,  <= 1.0 ,  := 0.03483533546735827;
	var t9 >= -1.0 ,  <= 1.0 ,  := 0.03108049841353322;
	var t10 >= -1.0 ,  <= 1.0 ,  := 0.027015115293406985;
	var x0 >= 0.0 ,  <= 0.0 ,  := 0.05;
	var x1 >= -0.05 ,  <= 0.05 ,  := 0.04975020826390129;
	var x2 >= -0.05 ,  <= 0.05 ,  := 0.04900332889206208;
	var x3 >= -0.05 ,  <= 0.05 ,  := 0.047766824456280305;
	var x4 >= -0.05 ,  <= 0.05 ,  := 0.04605304970014426;
	var x5 >= -0.05 ,  <= 0.05 ,  := 0.04387912809451864;
	var x6 >= -0.05 ,  <= 0.05 ,  := 0.041266780745483914;
	var x7 >= -0.05 ,  <= 0.05 ,  := 0.038242109364224425;
	var x8 >= -0.05 ,  <= 0.05 ,  := 0.03483533546735827;
	var x9 >= -0.05 ,  <= 0.05 ,  := 0.03108049841353322;
	var x10 >= 0.0 ,  <= 0.0 ,  := 0.027015115293406985;
	var u0;
	var u1;
	var u2;
	var u3;
	var u4;
	var u5;
	var u6;
	var u7;
	var u8;
	var u9;
	var u10;

minimize obj:
	0.05*u1 * u1 + 0.05*u0 * u0 + 17.5*(cos(t1)) + 17.5*(cos(t0)) + 0.05*u2 * u2 + 
	0.05*u1 * u1 + 17.5*(cos(t2)) + 17.5*(cos(t1)) + 0.05*u3 * u3 + 0.05*u2 * u2 + 
	17.5*(cos(t3)) + 17.5*(cos(t2)) + 0.05*u4 * u4 + 0.05*u3 * u3 + 17.5*(cos(t4)) 
	+ 17.5*(cos(t3)) + 0.05*u5 * u5 + 0.05*u4 * u4 + 17.5*(cos(t5)) + 
	17.5*(cos(t4)) + 0.05*u6 * u6 + 0.05*u5 * u5 + 17.5*(cos(t6)) + 17.5*(cos(t5)) 
	+ 0.05*u7 * u7 + 0.05*u6 * u6 + 17.5*(cos(t7)) + 17.5*(cos(t6)) + 0.05*u8 * u8 
	+ 0.05*u7 * u7 + 17.5*(cos(t8)) + 17.5*(cos(t7)) + 0.05*u9 * u9 + 0.05*u8 * u8 
	+ 17.5*(cos(t9)) + 17.5*(cos(t8)) + 0.05*u10 * u10 + 0.05*u9 * u9 + 
	17.5*(cos(t10)) + 17.5*(cos(t9));

subject to ex0:
	-0.05*(sin(t1)) - 0.05*(sin(t0)) + x1 - x0 = 0;
subject to et0:
	t1 - t0 - 0.05*u1 - 0.05*u0 = 0;
subject to ex1:
	-0.05*(sin(t2)) - 0.05*(sin(t1)) + x2 - x1 = 0;
subject to et1:
	t2 - t1 - 0.05*u2 - 0.05*u1 = 0;
subject to ex2:
	-0.05*(sin(t3)) - 0.05*(sin(t2)) + x3 - x2 = 0;
subject to et2:
	t3 - t2 - 0.05*u3 - 0.05*u2 = 0;
subject to ex3:
	-0.05*(sin(t4)) - 0.05*(sin(t3)) + x4 - x3 = 0;
subject to et3:
	t4 - t3 - 0.05*u4 - 0.05*u3 = 0;
subject to ex4:
	-0.05*(sin(t5)) - 0.05*(sin(t4)) + x5 - x4 = 0;
subject to et4:
	t5 - t4 - 0.05*u5 - 0.05*u4 = 0;
subject to ex5:
	-0.05*(sin(t6)) - 0.05*(sin(t5)) + x6 - x5 = 0;
subject to et5:
	t6 - t5 - 0.05*u6 - 0.05*u5 = 0;
subject to ex6:
	-0.05*(sin(t7)) - 0.05*(sin(t6)) + x7 - x6 = 0;
subject to et6:
	t7 - t6 - 0.05*u7 - 0.05*u6 = 0;
subject to ex7:
	-0.05*(sin(t8)) - 0.05*(sin(t7)) + x8 - x7 = 0;
subject to et7:
	t8 - t7 - 0.05*u8 - 0.05*u7 = 0;
subject to ex8:
	-0.05*(sin(t9)) - 0.05*(sin(t8)) + x9 - x8 = 0;
subject to et8:
	t9 - t8 - 0.05*u9 - 0.05*u8 = 0;
subject to ex9:
	-0.05*(sin(t10)) - 0.05*(sin(t9)) + x10 - x9 = 0;
subject to et9:
	t10 - t9 - 0.05*u10 - 0.05*u9 = 0;

solve;
	display t0;
	display t1;
	display t2;
	display t3;
	display t4;
	display t5;
	display t6;
	display t7;
	display t8;
	display t9;
	display t10;
	display x0;
	display x1;
	display x2;
	display x3;
	display x4;
	display x5;
	display x6;
	display x7;
	display x8;
	display x9;
	display x10;
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
	display u10;
display obj;
