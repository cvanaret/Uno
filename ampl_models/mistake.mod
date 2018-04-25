#
#**************************
# SET UP THE INITIAL DATA *
#**************************
#   Problem :
#   *********
#   A mistake in writing Hock and Schittkowski problem 108.
#   Source:
#   Ph. Toint.
#   classification QQR2-AY-9-13
#   SIF input: Ph. Toint, Apr 1990.
#   Number of variables
#   Parameters
#   Objective Function
#   Constraint function
#   Solution
	param n := 9;

	var x1 := 1.0;
	var x2 := 1.0;
	var x3 := 1.0;
	var x4 := 1.0;
	var x5 := 1.0;
	var x6 := 1.0;
	var x7 := 1.0;
	var x8 := 1.0;
	var x9 >= 0.0 ,  := 1.0;

minimize obj:
	 - 0.5*x1 * x4 + 0.5*x2 * x3 - 0.5*x3 * x9 + 0.5*x5 * x9 - 0.5*x5 * x8 + 0.5*x6 
	* x7;

subject to c1:
	0 >= x3 * x3 + x4 * x4 - 1.0;
subject to c2:
	0 >= x5 * x5 + x6 * x6 - 1.0;
subject to c3:
	0 >= x9 * x9 - 1.0;
subject to c4:
	0 >= x1 * x1 + (x2-x9) * (x2-x9) - 1.0;
subject to c5:
	0 >= (x1-x5) * (x1-x5) + (x2-x6) * (x2-x6) - 1.0;
subject to c6:
	0 >= (x1-x7) * (x1-x7) + (x2-x8) * (x2-x8) - 1.0;
subject to c7:
	0 >= (x3-x5) * (x3-x5) + (x4-x6) * (x4-x6) - 1.0;
subject to c8:
	0 >= (x3-x7) * (x3-x7) + (x4-x8) * (x4-x8) - 1.0;
subject to c9:
	0 >= x7 * x7 + x8 * x9 - 1.0;
subject to c10:
	0 <= x8 * x9;
subject to c11:
	0 <= x5 * x8 - x6 * x7;
subject to c12:
	0 <= x1 * x4 - x2 * x3;
subject to c13:
	0 >= -x5 * x9;

solve;
	display x1;
	display x2;
	display x3;
	display x4;
	display x5;
	display x6;
	display x7;
	display x8;
	display x9;
display obj;
