#
#**************************
# SET UP THE INITIAL DATA *
#**************************
#   Problem :
#   *********
#   Source: an expanded form of problem 99 in
#   W. Hock and K. Schittkowski,
#   "Test examples for nonlinear programming codes",
#   Lectures Notes in Economics and Mathematical Systems 187, Springer
#   Verlag, Heidelberg, 1981.
#   SIF input: Ph. Toint, April 1991.
#   classification OOR2-AN-31-21
#   Problem data
#   Constants
#   Increments
#   Solution
	param t1 := 0.0;
	param t2 := 25.0;
	param t3 := 50.0;
	param t4 := 100.0;
	param t5 := 150.0;
	param t6 := 200.0;
	param t7 := 290.0;
	param t8 := 380.0;
	param a1 := 0.0;
	param a2 := 50.0;
	param a3 := 50.0;
	param a4 := 75.0;
	param a5 := 75.0;
	param a6 := 75.0;
	param a7 := 100.0;
	param a8 := 100.0;
	param b := 32.0;
	param im1 := -1 + (8);
	param dtisq := ((380.0) - (290.0)) * ((380.0) - (290.0));
	param dt2 := (25.0) - (0.0);
	param dtsqd22 := 0.5 * (((25.0) - (0.0)) * ((25.0) - (0.0)));
	param dt3 := (50.0) - (25.0);
	param dtsqd23 := 0.5 * (((50.0) - (25.0)) * ((50.0) - (25.0)));
	param dt4 := (100.0) - (50.0);
	param dtsqd24 := 0.5 * (((100.0) - (50.0)) * ((100.0) - (50.0)));
	param dt5 := (150.0) - (100.0);
	param dtsqd25 := 0.5 * (((150.0) - (100.0)) * ((150.0) - (100.0)));
	param dt6 := (200.0) - (150.0);
	param dtsqd26 := 0.5 * (((200.0) - (150.0)) * ((200.0) - (150.0)));
	param dt7 := (290.0) - (200.0);
	param dtsqd27 := 0.5 * (((290.0) - (200.0)) * ((290.0) - (200.0)));
	param dt8 := (380.0) - (290.0);
	param dtsqd28 := 0.5 * (((380.0) - (290.0)) * ((380.0) - (290.0)));
	param rhs := ((290.0) - (200.0)) * (32.0);
	param w := (100.0) * (0.5 * (((380.0) - (290.0)) * ((380.0) - (290.0))));

	var x1 >= 0.0 ,  <= 1.58 ,  := 0.5;
	var r1 >= 0.0 ,  <= 0.0;
	var q1 >= 0.0 ,  <= 0.0;
	var s1 >= 0.0 ,  <= 0.0;
	var x2 >= 0.0 ,  <= 1.58 ,  := 0.5;
	var r2;
	var q2;
	var s2;
	var x3 >= 0.0 ,  <= 1.58 ,  := 0.5;
	var r3;
	var q3;
	var s3;
	var x4 >= 0.0 ,  <= 1.58 ,  := 0.5;
	var r4;
	var q4;
	var s4;
	var x5 >= 0.0 ,  <= 1.58 ,  := 0.5;
	var r5;
	var q5;
	var s5;
	var x6 >= 0.0 ,  <= 1.58 ,  := 0.5;
	var r6;
	var q6;
	var s6;
	var x7 >= 0.0 ,  <= 1.58 ,  := 0.5;
	var r7;
	var q7;
	var s7;
	var r8;
	var q8;
	var s8;

minimize obj:
	 - (r8)*(r8);

subject to rdef2:
	1250.0*(cos(x1)) - r2 + r1 = 0;
subject to qdef2:
	15625.0*(sin(x1)) - q2 + q1 + 25.0*s1 - 10000.0 = 0;
subject to sdef2:
	1250.0*(sin(x1)) - s2 + s1 - 800.0 = 0;
subject to rdef3:
	1250.0*(cos(x2)) - r3 + r2 = 0;
subject to qdef3:
	15625.0*(sin(x2)) - q3 + q2 + 25.0*s2 - 10000.0 = 0;
subject to sdef3:
	1250.0*(sin(x2)) - s3 + s2 - 800.0 = 0;
subject to rdef4:
	3750.0*(cos(x3)) - r4 + r3 = 0;
subject to qdef4:
	93750.0*(sin(x3)) - q4 + q3 + 50.0*s3 - 40000.0 = 0;
subject to sdef4:
	3750.0*(sin(x3)) - s4 + s3 - 1600.0 = 0;
subject to rdef5:
	3750.0*(cos(x4)) - r5 + r4 = 0;
subject to qdef5:
	93750.0*(sin(x4)) - q5 + q4 + 50.0*s4 - 40000.0 = 0;
subject to sdef5:
	3750.0*(sin(x4)) - s5 + s4 - 1600.0 = 0;
subject to rdef6:
	3750.0*(cos(x5)) - r6 + r5 = 0;
subject to qdef6:
	93750.0*(sin(x5)) - q6 + q5 + 50.0*s5 - 40000.0 = 0;
subject to sdef6:
	3750.0*(sin(x5)) - s6 + s5 - 1600.0 = 0;
subject to rdef7:
	9000.0*(cos(x6)) - r7 + r6 = 0;
subject to qdef7:
	405000.0*(sin(x6)) - q7 + q6 + 90.0*s6 - 129600.0 = 0;
subject to sdef7:
	9000.0*(sin(x6)) - s7 + s6 - 2880.0 = 0;
subject to rdef8:
	9000.0*(cos(x7)) - r8 + r7 = 0;
subject to qdef8:
	405000.0*(sin(x7)) - q8 + q7 + 90.0*s7 - 100000.0 = 0;
subject to sdef8:
	9000.0*(sin(x7)) - s8 + s7 - 1000.0 = 0;

solve;
	display x1;
	display r1;
	display q1;
	display s1;
	display x2;
	display r2;
	display q2;
	display s2;
	display x3;
	display r3;
	display q3;
	display s3;
	display x4;
	display r4;
	display q4;
	display s4;
	display x5;
	display r5;
	display q5;
	display s5;
	display x6;
	display r6;
	display q6;
	display s6;
	display x7;
	display r7;
	display q7;
	display s7;
	display r8;
	display q8;
	display s8;
display obj;
