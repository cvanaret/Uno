#
#**************************
# SET UP THE INITIAL DATA *
#**************************
#   Problem :
#   *********
#   Convex Square Root formulation of a nonconvex MINLP arising from
#   trim loss minimization in the paper industry. The problem is to 
#   produce a set of product paper rolls from raw paper rolls such that 
#   a cost function including the trim loss and the overall production cost
#   is minimized.
#
#   SMALLER PROBLEM
#   Source: I. Harjunkoski, T. Westerlund, R. P\"{o}rn and H. Skrifvars
#   "Different transformations for solving non--convex trim loss problems 
#   by MINLP", Report, Process Design Laboratory at Abo Akademi University,
#   Finland, Jan. 1996.
#   SIF input: S. Leyffer, November 1997
#   classification LOR2-AN-142-75
# PROBLEM SPECIFI CONSTANTS
# ... cost coefficients
# ... objective / cost function (1)
# ... constraint (2)
# ... constraint set (3)
# ... constraint set (4)
# ... constraint set (5)
# ... constraint set (6)
# ... constraint set (7)
# ... constraint set (8)
# ... constraint set (9)
# ... constraint set (12)
# ... constraint set (13)
# ... constraint set (14)
# ... constraint set (15)
# ... constraint set (16)
# ... nonlinear constraints (17) LINEAR PART
	param ii := 4;
	param jj := 4;
	param nmax := 5;
	param nord1 := 12.0;
	param nord2 := 6.0;
	param nord3 := 15.0;
	param nord4 := 9.0;
	param ll1 := 15.0;
	param ll2 := 12.0;
	param ll3 := 9.0;
	param ll4 := 6.0;
	param kk1 := 5.0;
	param kk2 := 5.0;
	param kk3 := 5.0;
	param kk4 := 4.0;
	param ten := 10.0;
	param bmax := 2000.0;
	param delta := 200.0;
	param b1 := 330.0;
	param b2 := 365.0;
	param b3 := 390.0;
	param b4 := 435.0;
	param count := precision(4, 2);
	param cc1 := (precision(1, 2)) / (10.0);
	param c1 := 1.0;
	param cc2 := (precision(2, 2)) / (10.0);
	param c2 := 1.0;
	param cc3 := (precision(3, 2)) / (10.0);
	param c3 := 1.0;
	param cc4 := (precision(4, 2)) / (10.0);
	param c4 := 1.0;
	param lljr := 6.0;
	param llj := round(6.0);
	param ind := precision(6, 2);
	param ccc1_1 := (precision(1, 2)) * (1.0);
	param ccc1_2 := (precision(2, 2)) * (1.0);
	param ccc1_3 := (precision(3, 2)) * (1.0);
	param ccc1_4 := (precision(4, 2)) * (1.0);
	param ccc1_5 := (precision(5, 2)) * (1.0);
	param ccc1_6 := (precision(6, 2)) * (1.0);
	param ccc1_7 := (precision(7, 2)) * (1.0);
	param ccc1_8 := (precision(8, 2)) * (1.0);
	param ccc1_9 := (precision(9, 2)) * (1.0);
	param ccc1_10 := (precision(10, 2)) * (1.0);
	param ccc1_11 := (precision(11, 2)) * (1.0);
	param ccc1_12 := (precision(12, 2)) * (1.0);
	param ccc1_13 := (precision(13, 2)) * (1.0);
	param ccc1_14 := (precision(14, 2)) * (1.0);
	param ccc1_15 := (precision(15, 2)) * (1.0);
	param ccc2_1 := (precision(1, 2)) * (1.0);
	param ccc2_2 := (precision(2, 2)) * (1.0);
	param ccc2_3 := (precision(3, 2)) * (1.0);
	param ccc2_4 := (precision(4, 2)) * (1.0);
	param ccc2_5 := (precision(5, 2)) * (1.0);
	param ccc2_6 := (precision(6, 2)) * (1.0);
	param ccc2_7 := (precision(7, 2)) * (1.0);
	param ccc2_8 := (precision(8, 2)) * (1.0);
	param ccc2_9 := (precision(9, 2)) * (1.0);
	param ccc2_10 := (precision(10, 2)) * (1.0);
	param ccc2_11 := (precision(11, 2)) * (1.0);
	param ccc2_12 := (precision(12, 2)) * (1.0);
	param ccc3_1 := (precision(1, 2)) * (1.0);
	param ccc3_2 := (precision(2, 2)) * (1.0);
	param ccc3_3 := (precision(3, 2)) * (1.0);
	param ccc3_4 := (precision(4, 2)) * (1.0);
	param ccc3_5 := (precision(5, 2)) * (1.0);
	param ccc3_6 := (precision(6, 2)) * (1.0);
	param ccc3_7 := (precision(7, 2)) * (1.0);
	param ccc3_8 := (precision(8, 2)) * (1.0);
	param ccc3_9 := (precision(9, 2)) * (1.0);
	param ccc4_1 := (precision(1, 2)) * (1.0);
	param ccc4_2 := (precision(2, 2)) * (1.0);
	param ccc4_3 := (precision(3, 2)) * (1.0);
	param ccc4_4 := (precision(4, 2)) * (1.0);
	param ccc4_5 := (precision(5, 2)) * (1.0);
	param ccc4_6 := (precision(6, 2)) * (1.0);
	param kkir := 4.0;
	param kki := round(4.0);
	param temp := 4;
	param lr := precision(6, 2);
	param kr := precision(4, 2);
	param jp1 := 1 + (3);
	param lljp1r := 6.0;
	param lljp1 := round(6.0);
	param rhs := ((0.0) - (9.0)) - (4);

	var y1 >= 0.0 ,  <= 1.0;
	var y2 >= 0.0 ,  <= 1.0;
	var y3 >= 0.0 ,  <= 1.0;
	var y4 >= 0.0 ,  <= 1.0;
	var beta1_1 >= 0.0 ,  <= 1.0;
	var beta1_2 >= 0.0 ,  <= 1.0;
	var beta1_3 >= 0.0 ,  <= 1.0;
	var beta1_4 >= 0.0 ,  <= 1.0;
	var beta1_5 >= 0.0 ,  <= 1.0;
	var beta1_6 >= 0.0 ,  <= 1.0;
	var beta1_7 >= 0.0 ,  <= 1.0;
	var beta1_8 >= 0.0 ,  <= 1.0;
	var beta1_9 >= 0.0 ,  <= 1.0;
	var beta1_10 >= 0.0 ,  <= 1.0;
	var beta1_11 >= 0.0 ,  <= 1.0;
	var beta1_12 >= 0.0 ,  <= 1.0;
	var beta1_13 >= 0.0 ,  <= 1.0;
	var beta1_14 >= 0.0 ,  <= 1.0;
	var beta1_15 >= 0.0 ,  <= 1.0;
	var beta2_1 >= 0.0 ,  <= 1.0;
	var beta2_2 >= 0.0 ,  <= 1.0;
	var beta2_3 >= 0.0 ,  <= 1.0;
	var beta2_4 >= 0.0 ,  <= 1.0;
	var beta2_5 >= 0.0 ,  <= 1.0;
	var beta2_6 >= 0.0 ,  <= 1.0;
	var beta2_7 >= 0.0 ,  <= 1.0;
	var beta2_8 >= 0.0 ,  <= 1.0;
	var beta2_9 >= 0.0 ,  <= 1.0;
	var beta2_10 >= 0.0 ,  <= 1.0;
	var beta2_11 >= 0.0 ,  <= 1.0;
	var beta2_12 >= 0.0 ,  <= 1.0;
	var beta3_1 >= 0.0 ,  <= 1.0;
	var beta3_2 >= 0.0 ,  <= 1.0;
	var beta3_3 >= 0.0 ,  <= 1.0;
	var beta3_4 >= 0.0 ,  <= 1.0;
	var beta3_5 >= 0.0 ,  <= 1.0;
	var beta3_6 >= 0.0 ,  <= 1.0;
	var beta3_7 >= 0.0 ,  <= 1.0;
	var beta3_8 >= 0.0 ,  <= 1.0;
	var beta3_9 >= 0.0 ,  <= 1.0;
	var beta4_1 >= 0.0 ,  <= 1.0;
	var beta4_2 >= 0.0 ,  <= 1.0;
	var beta4_3 >= 0.0 ,  <= 1.0;
	var beta4_4 >= 0.0 ,  <= 1.0;
	var beta4_5 >= 0.0 ,  <= 1.0;
	var beta4_6 >= 0.0 ,  <= 1.0;
	var z1_1_1 >= 0.0 ,  <= 1.0;
	var z1_1_2 >= 0.0 ,  <= 1.0;
	var z1_1_3 >= 0.0 ,  <= 1.0;
	var z1_1_4 >= 0.0 ,  <= 1.0;
	var z1_1_5 >= 0.0 ,  <= 1.0;
	var z1_2_1 >= 0.0 ,  <= 1.0;
	var z1_2_2 >= 0.0 ,  <= 1.0;
	var z1_2_3 >= 0.0 ,  <= 1.0;
	var z1_2_4 >= 0.0 ,  <= 1.0;
	var z1_2_5 >= 0.0 ,  <= 1.0;
	var z1_3_1 >= 0.0 ,  <= 1.0;
	var z1_3_2 >= 0.0 ,  <= 1.0;
	var z1_3_3 >= 0.0 ,  <= 1.0;
	var z1_3_4 >= 0.0 ,  <= 1.0;
	var z1_3_5 >= 0.0 ,  <= 1.0;
	var z1_4_1 >= 0.0 ,  <= 1.0;
	var z1_4_2 >= 0.0 ,  <= 1.0;
	var z1_4_3 >= 0.0 ,  <= 1.0;
	var z1_4_4 >= 0.0 ,  <= 1.0;
	var z1_4_5 >= 0.0 ,  <= 1.0;
	var z2_1_1 >= 0.0 ,  <= 1.0;
	var z2_1_2 >= 0.0 ,  <= 1.0;
	var z2_1_3 >= 0.0 ,  <= 1.0;
	var z2_1_4 >= 0.0 ,  <= 1.0;
	var z2_1_5 >= 0.0 ,  <= 1.0;
	var z2_2_1 >= 0.0 ,  <= 1.0;
	var z2_2_2 >= 0.0 ,  <= 1.0;
	var z2_2_3 >= 0.0 ,  <= 1.0;
	var z2_2_4 >= 0.0 ,  <= 1.0;
	var z2_2_5 >= 0.0 ,  <= 1.0;
	var z2_3_1 >= 0.0 ,  <= 1.0;
	var z2_3_2 >= 0.0 ,  <= 1.0;
	var z2_3_3 >= 0.0 ,  <= 1.0;
	var z2_3_4 >= 0.0 ,  <= 1.0;
	var z2_3_5 >= 0.0 ,  <= 1.0;
	var z2_4_1 >= 0.0 ,  <= 1.0;
	var z2_4_2 >= 0.0 ,  <= 1.0;
	var z2_4_3 >= 0.0 ,  <= 1.0;
	var z2_4_4 >= 0.0 ,  <= 1.0;
	var z2_4_5 >= 0.0 ,  <= 1.0;
	var z3_1_1 >= 0.0 ,  <= 1.0;
	var z3_1_2 >= 0.0 ,  <= 1.0;
	var z3_1_3 >= 0.0 ,  <= 1.0;
	var z3_1_4 >= 0.0 ,  <= 1.0;
	var z3_1_5 >= 0.0 ,  <= 1.0;
	var z3_2_1 >= 0.0 ,  <= 1.0;
	var z3_2_2 >= 0.0 ,  <= 1.0;
	var z3_2_3 >= 0.0 ,  <= 1.0;
	var z3_2_4 >= 0.0 ,  <= 1.0;
	var z3_2_5 >= 0.0 ,  <= 1.0;
	var z3_3_1 >= 0.0 ,  <= 1.0;
	var z3_3_2 >= 0.0 ,  <= 1.0;
	var z3_3_3 >= 0.0 ,  <= 1.0;
	var z3_3_4 >= 0.0 ,  <= 1.0;
	var z3_3_5 >= 0.0 ,  <= 1.0;
	var z3_4_1 >= 0.0 ,  <= 1.0;
	var z3_4_2 >= 0.0 ,  <= 1.0;
	var z3_4_3 >= 0.0 ,  <= 1.0;
	var z3_4_4 >= 0.0 ,  <= 1.0;
	var z3_4_5 >= 0.0 ,  <= 1.0;
	var z4_1_1 >= 0.0 ,  <= 1.0;
	var z4_1_2 >= 0.0 ,  <= 1.0;
	var z4_1_3 >= 0.0 ,  <= 1.0;
	var z4_1_4 >= 0.0 ,  <= 1.0;
	var z4_2_1 >= 0.0 ,  <= 1.0;
	var z4_2_2 >= 0.0 ,  <= 1.0;
	var z4_2_3 >= 0.0 ,  <= 1.0;
	var z4_2_4 >= 0.0 ,  <= 1.0;
	var z4_3_1 >= 0.0 ,  <= 1.0;
	var z4_3_2 >= 0.0 ,  <= 1.0;
	var z4_3_3 >= 0.0 ,  <= 1.0;
	var z4_3_4 >= 0.0 ,  <= 1.0;
	var z4_4_1 >= 0.0 ,  <= 1.0;
	var z4_4_2 >= 0.0 ,  <= 1.0;
	var z4_4_3 >= 0.0 ,  <= 1.0;
	var z4_4_4 >= 0.0 ,  <= 1.0;
	var m1 >= 0.0 ,  := 1.0;
	var m2 >= 0.0 ,  := 1.0;
	var m3 >= 0.0 ,  := 1.0;
	var m4 >= 0.0 ,  := 1.0;
	var n1_1 >= 0.0 ,  := 1.0;
	var n1_2 >= 0.0 ,  := 1.0;
	var n1_3 >= 0.0 ,  := 1.0;
	var n1_4 >= 0.0 ,  := 1.0;
	var n2_1 >= 0.0 ,  := 1.0;
	var n2_2 >= 0.0 ,  := 1.0;
	var n2_3 >= 0.0 ,  := 1.0;
	var n2_4 >= 0.0 ,  := 1.0;
	var n3_1 >= 0.0 ,  := 1.0;
	var n3_2 >= 0.0 ,  := 1.0;
	var n3_3 >= 0.0 ,  := 1.0;
	var n3_4 >= 0.0 ,  := 1.0;
	var n4_1 >= 0.0 ,  := 1.0;
	var n4_2 >= 0.0 ,  := 1.0;
	var n4_3 >= 0.0 ,  := 1.0;
	var n4_4 >= 0.0 ,  := 1.0;

minimize obj:
	0.1*y1 + 0.2*y2 + 0.3*y3 + 0.4*y4 + beta1_1 + 2.0*beta1_2 + 3.0*beta1_3 + 
	4.0*beta1_4 + 5.0*beta1_5 + 6.0*beta1_6 + 7.0*beta1_7 + 8.0*beta1_8 + 
	9.0*beta1_9 + 10.0*beta1_10 + 11.0*beta1_11 + 12.0*beta1_12 + 13.0*beta1_13 + 
	14.0*beta1_14 + 15.0*beta1_15 + beta2_1 + 2.0*beta2_2 + 3.0*beta2_3 + 
	4.0*beta2_4 + 5.0*beta2_5 + 6.0*beta2_6 + 7.0*beta2_7 + 8.0*beta2_8 + 
	9.0*beta2_9 + 10.0*beta2_10 + 11.0*beta2_11 + 12.0*beta2_12 + beta3_1 + 
	2.0*beta3_2 + 3.0*beta3_3 + 4.0*beta3_4 + 5.0*beta3_5 + 6.0*beta3_6 + 
	7.0*beta3_7 + 8.0*beta3_8 + 9.0*beta3_9 + beta4_1 + 2.0*beta4_2 + 3.0*beta4_3 + 
	4.0*beta4_4 + 5.0*beta4_5 + 6.0*beta4_6;

subject to lin2:
	0 >= -beta1_1 - 2.0*beta1_2 - 3.0*beta1_3 - 4.0*beta1_4 - 5.0*beta1_5 - 
	6.0*beta1_6 - 7.0*beta1_7 - 8.0*beta1_8 - 9.0*beta1_9 - 10.0*beta1_10 - 
	11.0*beta1_11 - 12.0*beta1_12 - 13.0*beta1_13 - 14.0*beta1_14 - 15.0*beta1_15 - 
	beta2_1 - 2.0*beta2_2 - 3.0*beta2_3 - 4.0*beta2_4 - 5.0*beta2_5 - 6.0*beta2_6 - 
	7.0*beta2_7 - 8.0*beta2_8 - 9.0*beta2_9 - 10.0*beta2_10 - 11.0*beta2_11 - 
	12.0*beta2_12 - beta3_1 - 2.0*beta3_2 - 3.0*beta3_3 - 4.0*beta3_4 - 5.0*beta3_5 
	- 6.0*beta3_6 - 7.0*beta3_7 - 8.0*beta3_8 - 9.0*beta3_9 - beta4_1 - 2.0*beta4_2 
	- 3.0*beta4_3 - 4.0*beta4_4 - 5.0*beta4_5 - 6.0*beta4_6 + 9.0;
subject to lin31:
	0 >= z1_1_1 + 2.0*z1_1_2 + 3.0*z1_1_3 + 4.0*z1_1_4 + 5.0*z1_1_5 - 5.0*y1 + 
	z2_1_1 + 2.0*z2_1_2 + 3.0*z2_1_3 + 4.0*z2_1_4 + 5.0*z2_1_5 - 5.0*y1 + z3_1_1 + 
	2.0*z3_1_2 + 3.0*z3_1_3 + 4.0*z3_1_4 + 5.0*z3_1_5 - 5.0*y1 + z4_1_1 + 
	2.0*z4_1_2 + 3.0*z4_1_3 + 4.0*z4_1_4 - 4.0*y1;
subject to lin32:
	0 >= z1_2_1 + 2.0*z1_2_2 + 3.0*z1_2_3 + 4.0*z1_2_4 + 5.0*z1_2_5 - 5.0*y2 + 
	z2_2_1 + 2.0*z2_2_2 + 3.0*z2_2_3 + 4.0*z2_2_4 + 5.0*z2_2_5 - 5.0*y2 + z3_2_1 + 
	2.0*z3_2_2 + 3.0*z3_2_3 + 4.0*z3_2_4 + 5.0*z3_2_5 - 5.0*y2 + z4_2_1 + 
	2.0*z4_2_2 + 3.0*z4_2_3 + 4.0*z4_2_4 - 4.0*y2;
subject to lin33:
	0 >= z1_3_1 + 2.0*z1_3_2 + 3.0*z1_3_3 + 4.0*z1_3_4 + 5.0*z1_3_5 - 5.0*y3 + 
	z2_3_1 + 2.0*z2_3_2 + 3.0*z2_3_3 + 4.0*z2_3_4 + 5.0*z2_3_5 - 5.0*y3 + z3_3_1 + 
	2.0*z3_3_2 + 3.0*z3_3_3 + 4.0*z3_3_4 + 5.0*z3_3_5 - 5.0*y3 + z4_3_1 + 
	2.0*z4_3_2 + 3.0*z4_3_3 + 4.0*z4_3_4 - 4.0*y3;
subject to lin34:
	0 >= z1_4_1 + 2.0*z1_4_2 + 3.0*z1_4_3 + 4.0*z1_4_4 + 5.0*z1_4_5 - 5.0*y4 + 
	z2_4_1 + 2.0*z2_4_2 + 3.0*z2_4_3 + 4.0*z2_4_4 + 5.0*z2_4_5 - 5.0*y4 + z3_4_1 + 
	2.0*z3_4_2 + 3.0*z3_4_3 + 4.0*z3_4_4 + 5.0*z3_4_5 - 5.0*y4 + z4_4_1 + 
	2.0*z4_4_2 + 3.0*z4_4_3 + 4.0*z4_4_4 - 4.0*y4;
subject to lin41:
	0 >= 330.0*z1_1_1 + 660.0*z1_1_2 + 990.0*z1_1_3 + 1320.0*z1_1_4 + 1650.0*z1_1_5 
	+ 365.0*z2_1_1 + 730.0*z2_1_2 + 1095.0*z2_1_3 + 1460.0*z2_1_4 + 1825.0*z2_1_5 + 
	390.0*z3_1_1 + 780.0*z3_1_2 + 1170.0*z3_1_3 + 1560.0*z3_1_4 + 1950.0*z3_1_5 + 
	435.0*z4_1_1 + 870.0*z4_1_2 + 1305.0*z4_1_3 + 1740.0*z4_1_4 - 2000.0*y1;
subject to lin42:
	0 >= 330.0*z1_2_1 + 660.0*z1_2_2 + 990.0*z1_2_3 + 1320.0*z1_2_4 + 1650.0*z1_2_5 
	+ 365.0*z2_2_1 + 730.0*z2_2_2 + 1095.0*z2_2_3 + 1460.0*z2_2_4 + 1825.0*z2_2_5 + 
	390.0*z3_2_1 + 780.0*z3_2_2 + 1170.0*z3_2_3 + 1560.0*z3_2_4 + 1950.0*z3_2_5 + 
	435.0*z4_2_1 + 870.0*z4_2_2 + 1305.0*z4_2_3 + 1740.0*z4_2_4 - 2000.0*y2;
subject to lin43:
	0 >= 330.0*z1_3_1 + 660.0*z1_3_2 + 990.0*z1_3_3 + 1320.0*z1_3_4 + 1650.0*z1_3_5 
	+ 365.0*z2_3_1 + 730.0*z2_3_2 + 1095.0*z2_3_3 + 1460.0*z2_3_4 + 1825.0*z2_3_5 + 
	390.0*z3_3_1 + 780.0*z3_3_2 + 1170.0*z3_3_3 + 1560.0*z3_3_4 + 1950.0*z3_3_5 + 
	435.0*z4_3_1 + 870.0*z4_3_2 + 1305.0*z4_3_3 + 1740.0*z4_3_4 - 2000.0*y3;
subject to lin44:
	0 >= 330.0*z1_4_1 + 660.0*z1_4_2 + 990.0*z1_4_3 + 1320.0*z1_4_4 + 1650.0*z1_4_5 
	+ 365.0*z2_4_1 + 730.0*z2_4_2 + 1095.0*z2_4_3 + 1460.0*z2_4_4 + 1825.0*z2_4_5 + 
	390.0*z3_4_1 + 780.0*z3_4_2 + 1170.0*z3_4_3 + 1560.0*z3_4_4 + 1950.0*z3_4_5 + 
	435.0*z4_4_1 + 870.0*z4_4_2 + 1305.0*z4_4_3 + 1740.0*z4_4_4 - 2000.0*y4;
subject to lin51:
	0 >= -330.0*z1_1_1 - 660.0*z1_1_2 - 990.0*z1_1_3 - 1320.0*z1_1_4 - 
	1650.0*z1_1_5 - 365.0*z2_1_1 - 730.0*z2_1_2 - 1095.0*z2_1_3 - 1460.0*z2_1_4 - 
	1825.0*z2_1_5 - 390.0*z3_1_1 - 780.0*z3_1_2 - 1170.0*z3_1_3 - 1560.0*z3_1_4 - 
	1950.0*z3_1_5 - 435.0*z4_1_1 - 870.0*z4_1_2 - 1305.0*z4_1_3 - 1740.0*z4_1_4 + 
	1800.0*y1;
subject to lin52:
	0 >= -330.0*z1_2_1 - 660.0*z1_2_2 - 990.0*z1_2_3 - 1320.0*z1_2_4 - 
	1650.0*z1_2_5 - 365.0*z2_2_1 - 730.0*z2_2_2 - 1095.0*z2_2_3 - 1460.0*z2_2_4 - 
	1825.0*z2_2_5 - 390.0*z3_2_1 - 780.0*z3_2_2 - 1170.0*z3_2_3 - 1560.0*z3_2_4 - 
	1950.0*z3_2_5 - 435.0*z4_2_1 - 870.0*z4_2_2 - 1305.0*z4_2_3 - 1740.0*z4_2_4 + 
	1800.0*y2;
subject to lin53:
	0 >= -330.0*z1_3_1 - 660.0*z1_3_2 - 990.0*z1_3_3 - 1320.0*z1_3_4 - 
	1650.0*z1_3_5 - 365.0*z2_3_1 - 730.0*z2_3_2 - 1095.0*z2_3_3 - 1460.0*z2_3_4 - 
	1825.0*z2_3_5 - 390.0*z3_3_1 - 780.0*z3_3_2 - 1170.0*z3_3_3 - 1560.0*z3_3_4 - 
	1950.0*z3_3_5 - 435.0*z4_3_1 - 870.0*z4_3_2 - 1305.0*z4_3_3 - 1740.0*z4_3_4 + 
	1800.0*y3;
subject to lin54:
	0 >= -330.0*z1_4_1 - 660.0*z1_4_2 - 990.0*z1_4_3 - 1320.0*z1_4_4 - 
	1650.0*z1_4_5 - 365.0*z2_4_1 - 730.0*z2_4_2 - 1095.0*z2_4_3 - 1460.0*z2_4_4 - 
	1825.0*z2_4_5 - 390.0*z3_4_1 - 780.0*z3_4_2 - 1170.0*z3_4_3 - 1560.0*z3_4_4 - 
	1950.0*z3_4_5 - 435.0*z4_4_1 - 870.0*z4_4_2 - 1305.0*z4_4_3 - 1740.0*z4_4_4 + 
	1800.0*y4;
subject to lin61:
	0 >= -beta1_1 - 2.0*beta1_2 - 3.0*beta1_3 - 4.0*beta1_4 - 5.0*beta1_5 - 
	6.0*beta1_6 - 7.0*beta1_7 - 8.0*beta1_8 - 9.0*beta1_9 - 10.0*beta1_10 - 
	11.0*beta1_11 - 12.0*beta1_12 - 13.0*beta1_13 - 14.0*beta1_14 - 15.0*beta1_15 + 
	y1;
subject to lin62:
	0 >= -beta2_1 - 2.0*beta2_2 - 3.0*beta2_3 - 4.0*beta2_4 - 5.0*beta2_5 - 
	6.0*beta2_6 - 7.0*beta2_7 - 8.0*beta2_8 - 9.0*beta2_9 - 10.0*beta2_10 - 
	11.0*beta2_11 - 12.0*beta2_12 + y2;
subject to lin63:
	0 >= -beta3_1 - 2.0*beta3_2 - 3.0*beta3_3 - 4.0*beta3_4 - 5.0*beta3_5 - 
	6.0*beta3_6 - 7.0*beta3_7 - 8.0*beta3_8 - 9.0*beta3_9 + y3;
subject to lin64:
	0 >= -beta4_1 - 2.0*beta4_2 - 3.0*beta4_3 - 4.0*beta4_4 - 5.0*beta4_5 - 
	6.0*beta4_6 + y4;
subject to lin71:
	0 >= beta1_1 + 2.0*beta1_2 + 3.0*beta1_3 + 4.0*beta1_4 + 5.0*beta1_5 + 
	6.0*beta1_6 + 7.0*beta1_7 + 8.0*beta1_8 + 9.0*beta1_9 + 10.0*beta1_10 + 
	11.0*beta1_11 + 12.0*beta1_12 + 13.0*beta1_13 + 14.0*beta1_14 + 15.0*beta1_15 - 
	15.0*y1;
subject to lin72:
	0 >= beta2_1 + 2.0*beta2_2 + 3.0*beta2_3 + 4.0*beta2_4 + 5.0*beta2_5 + 
	6.0*beta2_6 + 7.0*beta2_7 + 8.0*beta2_8 + 9.0*beta2_9 + 10.0*beta2_10 + 
	11.0*beta2_11 + 12.0*beta2_12 - 12.0*y2;
subject to lin73:
	0 >= beta3_1 + 2.0*beta3_2 + 3.0*beta3_3 + 4.0*beta3_4 + 5.0*beta3_5 + 
	6.0*beta3_6 + 7.0*beta3_7 + 8.0*beta3_8 + 9.0*beta3_9 - 9.0*y3;
subject to lin74:
	0 >= beta4_1 + 2.0*beta4_2 + 3.0*beta4_3 + 4.0*beta4_4 + 5.0*beta4_5 + 
	6.0*beta4_6 - 6.0*y4;
subject to lin81:
	-3.0*beta1_1 - 8.0*beta1_2 - 15.0*beta1_3 - 24.0*beta1_4 - 35.0*beta1_5 - 
	48.0*beta1_6 - 63.0*beta1_7 - 80.0*beta1_8 - 99.0*beta1_9 - 120.0*beta1_10 - 
	143.0*beta1_11 - 168.0*beta1_12 - 195.0*beta1_13 - 224.0*beta1_14 - 
	255.0*beta1_15 + m1 - 1.0 = 0;
subject to lin82:
	-3.0*beta2_1 - 8.0*beta2_2 - 15.0*beta2_3 - 24.0*beta2_4 - 35.0*beta2_5 - 
	48.0*beta2_6 - 63.0*beta2_7 - 80.0*beta2_8 - 99.0*beta2_9 - 120.0*beta2_10 - 
	143.0*beta2_11 - 168.0*beta2_12 + m2 - 1.0 = 0;
subject to lin83:
	-3.0*beta3_1 - 8.0*beta3_2 - 15.0*beta3_3 - 24.0*beta3_4 - 35.0*beta3_5 - 
	48.0*beta3_6 - 63.0*beta3_7 - 80.0*beta3_8 - 99.0*beta3_9 + m3 - 1.0 = 0;
subject to lin84:
	-3.0*beta4_1 - 8.0*beta4_2 - 15.0*beta4_3 - 24.0*beta4_4 - 35.0*beta4_5 - 
	48.0*beta4_6 + m4 - 1.0 = 0;
subject to l91_1:
	-3.0*z1_1_1 - 8.0*z1_1_2 - 15.0*z1_1_3 - 24.0*z1_1_4 - 35.0*z1_1_5 + n1_1 - 1.0 
	= 0;
subject to l91_2:
	-3.0*z1_2_1 - 8.0*z1_2_2 - 15.0*z1_2_3 - 24.0*z1_2_4 - 35.0*z1_2_5 + n1_2 - 1.0 
	= 0;
subject to l91_3:
	-3.0*z1_3_1 - 8.0*z1_3_2 - 15.0*z1_3_3 - 24.0*z1_3_4 - 35.0*z1_3_5 + n1_3 - 1.0 
	= 0;
subject to l91_4:
	-3.0*z1_4_1 - 8.0*z1_4_2 - 15.0*z1_4_3 - 24.0*z1_4_4 - 35.0*z1_4_5 + n1_4 - 1.0 
	= 0;
subject to l92_1:
	-3.0*z2_1_1 - 8.0*z2_1_2 - 15.0*z2_1_3 - 24.0*z2_1_4 - 35.0*z2_1_5 + n2_1 - 1.0 
	= 0;
subject to l92_2:
	-3.0*z2_2_1 - 8.0*z2_2_2 - 15.0*z2_2_3 - 24.0*z2_2_4 - 35.0*z2_2_5 + n2_2 - 1.0 
	= 0;
subject to l92_3:
	-3.0*z2_3_1 - 8.0*z2_3_2 - 15.0*z2_3_3 - 24.0*z2_3_4 - 35.0*z2_3_5 + n2_3 - 1.0 
	= 0;
subject to l92_4:
	-3.0*z2_4_1 - 8.0*z2_4_2 - 15.0*z2_4_3 - 24.0*z2_4_4 - 35.0*z2_4_5 + n2_4 - 1.0 
	= 0;
subject to l93_1:
	-3.0*z3_1_1 - 8.0*z3_1_2 - 15.0*z3_1_3 - 24.0*z3_1_4 - 35.0*z3_1_5 + n3_1 - 1.0 
	= 0;
subject to l93_2:
	-3.0*z3_2_1 - 8.0*z3_2_2 - 15.0*z3_2_3 - 24.0*z3_2_4 - 35.0*z3_2_5 + n3_2 - 1.0 
	= 0;
subject to l93_3:
	-3.0*z3_3_1 - 8.0*z3_3_2 - 15.0*z3_3_3 - 24.0*z3_3_4 - 35.0*z3_3_5 + n3_3 - 1.0 
	= 0;
subject to l93_4:
	-3.0*z3_4_1 - 8.0*z3_4_2 - 15.0*z3_4_3 - 24.0*z3_4_4 - 35.0*z3_4_5 + n3_4 - 1.0 
	= 0;
subject to l94_1:
	-3.0*z4_1_1 - 8.0*z4_1_2 - 15.0*z4_1_3 - 24.0*z4_1_4 + n4_1 - 1.0 = 0;
subject to l94_2:
	-3.0*z4_2_1 - 8.0*z4_2_2 - 15.0*z4_2_3 - 24.0*z4_2_4 + n4_2 - 1.0 = 0;
subject to l94_3:
	-3.0*z4_3_1 - 8.0*z4_3_2 - 15.0*z4_3_3 - 24.0*z4_3_4 + n4_3 - 1.0 = 0;
subject to l94_4:
	-3.0*z4_4_1 - 8.0*z4_4_2 - 15.0*z4_4_3 - 24.0*z4_4_4 + n4_4 - 1.0 = 0;
subject to l121_1:
	0 >= z1_1_1 + z1_1_2 + z1_1_3 + z1_1_4 + z1_1_5 - 1.0;
subject to l121_2:
	0 >= z1_2_1 + z1_2_2 + z1_2_3 + z1_2_4 + z1_2_5 - 1.0;
subject to l121_3:
	0 >= z1_3_1 + z1_3_2 + z1_3_3 + z1_3_4 + z1_3_5 - 1.0;
subject to l121_4:
	0 >= z1_4_1 + z1_4_2 + z1_4_3 + z1_4_4 + z1_4_5 - 1.0;
subject to l122_1:
	0 >= z2_1_1 + z2_1_2 + z2_1_3 + z2_1_4 + z2_1_5 - 1.0;
subject to l122_2:
	0 >= z2_2_1 + z2_2_2 + z2_2_3 + z2_2_4 + z2_2_5 - 1.0;
subject to l122_3:
	0 >= z2_3_1 + z2_3_2 + z2_3_3 + z2_3_4 + z2_3_5 - 1.0;
subject to l122_4:
	0 >= z2_4_1 + z2_4_2 + z2_4_3 + z2_4_4 + z2_4_5 - 1.0;
subject to l123_1:
	0 >= z3_1_1 + z3_1_2 + z3_1_3 + z3_1_4 + z3_1_5 - 1.0;
subject to l123_2:
	0 >= z3_2_1 + z3_2_2 + z3_2_3 + z3_2_4 + z3_2_5 - 1.0;
subject to l123_3:
	0 >= z3_3_1 + z3_3_2 + z3_3_3 + z3_3_4 + z3_3_5 - 1.0;
subject to l123_4:
	0 >= z3_4_1 + z3_4_2 + z3_4_3 + z3_4_4 + z3_4_5 - 1.0;
subject to l124_1:
	0 >= z4_1_1 + z4_1_2 + z4_1_3 + z4_1_4 - 1.0;
subject to l124_2:
	0 >= z4_2_1 + z4_2_2 + z4_2_3 + z4_2_4 - 1.0;
subject to l124_3:
	0 >= z4_3_1 + z4_3_2 + z4_3_3 + z4_3_4 - 1.0;
subject to l124_4:
	0 >= z4_4_1 + z4_4_2 + z4_4_3 + z4_4_4 - 1.0;
subject to lin131:
	0 >= beta1_1 + beta1_2 + beta1_3 + beta1_4 + beta1_5 + beta1_6 + beta1_7 + 
	beta1_8 + beta1_9 + beta1_10 + beta1_11 + beta1_12 + beta1_13 + beta1_14 + 
	beta1_15 - 1.0;
subject to lin132:
	0 >= beta2_1 + beta2_2 + beta2_3 + beta2_4 + beta2_5 + beta2_6 + beta2_7 + 
	beta2_8 + beta2_9 + beta2_10 + beta2_11 + beta2_12 - 1.0;
subject to lin133:
	0 >= beta3_1 + beta3_2 + beta3_3 + beta3_4 + beta3_5 + beta3_6 + beta3_7 + 
	beta3_8 + beta3_9 - 1.0;
subject to lin134:
	0 >= beta4_1 + beta4_2 + beta4_3 + beta4_4 + beta4_5 + beta4_6 - 1.0;
subject to lin141:
	0 >= beta1_1 + beta1_2 + beta1_3 + beta1_4 + beta1_5 + beta1_6 + beta1_7 + 
	beta1_8 + beta1_9 + beta1_10 + beta1_11 + beta1_12 - beta1_1 - beta1_2 - 
	beta1_3 - beta1_4 - beta1_5 - beta1_6 - beta1_7 - beta1_8 - beta1_9 - beta1_10 
	- beta1_11 - beta1_12 - beta1_13 - beta1_14 - beta1_15;
subject to lin142:
	0 >= beta2_1 + beta2_2 + beta2_3 + beta2_4 + beta2_5 + beta2_6 + beta2_7 + 
	beta2_8 + beta2_9 - beta2_1 - beta2_2 - beta2_3 - beta2_4 - beta2_5 - beta2_6 - 
	beta2_7 - beta2_8 - beta2_9 - beta2_10 - beta2_11 - beta2_12;
subject to lin143:
	0 >= beta3_1 + beta3_2 + beta3_3 + beta3_4 + beta3_5 + beta3_6 - beta3_1 - 
	beta3_2 - beta3_3 - beta3_4 - beta3_5 - beta3_6 - beta3_7 - beta3_8 - beta3_9;
subject to lin151:
	0 >= y2 - y1;
subject to lin152:
	0 >= y3 - y2;
subject to lin153:
	0 >= y4 - y3;
subject to lin161:
	0 >= beta1_1 + 2.0*beta1_2 + 3.0*beta1_3 + 4.0*beta1_4 + 5.0*beta1_5 + 
	6.0*beta1_6 + 7.0*beta1_7 + 8.0*beta1_8 + 9.0*beta1_9 + 10.0*beta1_10 + 
	11.0*beta1_11 + 12.0*beta1_12 + 13.0*beta1_13 + 14.0*beta1_14 + 15.0*beta1_15 - 
	15.0;
subject to lin162:
	0 >= beta2_1 + 2.0*beta2_2 + 3.0*beta2_3 + 4.0*beta2_4 + 5.0*beta2_5 + 
	6.0*beta2_6 + 7.0*beta2_7 + 8.0*beta2_8 + 9.0*beta2_9 + 10.0*beta2_10 + 
	11.0*beta2_11 + 12.0*beta2_12 - 12.0;
subject to lin163:
	0 >= beta3_1 + 2.0*beta3_2 + 3.0*beta3_3 + 4.0*beta3_4 + 5.0*beta3_5 + 
	6.0*beta3_6 + 7.0*beta3_7 + 8.0*beta3_8 + 9.0*beta3_9 - 9.0;
subject to lin164:
	0 >= beta4_1 + 2.0*beta4_2 + 3.0*beta4_3 + 4.0*beta4_4 + 5.0*beta4_5 + 
	6.0*beta4_6 - 6.0;
subject to nln171:
	0 >= -(sqrt((m1*n1_1))) - (sqrt((m2*n1_2))) - (sqrt((m3*n1_3))) - 
	(sqrt((m4*n1_4))) + beta1_1 + 2.0*beta1_2 + 3.0*beta1_3 + 4.0*beta1_4 + 
	5.0*beta1_5 + 6.0*beta1_6 + 7.0*beta1_7 + 8.0*beta1_8 + 9.0*beta1_9 + 
	10.0*beta1_10 + 11.0*beta1_11 + 12.0*beta1_12 + 13.0*beta1_13 + 14.0*beta1_14 + 
	15.0*beta1_15 + z1_1_1 + 2.0*z1_1_2 + 3.0*z1_1_3 + 4.0*z1_1_4 + 5.0*z1_1_5 + 
	beta2_1 + 2.0*beta2_2 + 3.0*beta2_3 + 4.0*beta2_4 + 5.0*beta2_5 + 6.0*beta2_6 + 
	7.0*beta2_7 + 8.0*beta2_8 + 9.0*beta2_9 + 10.0*beta2_10 + 11.0*beta2_11 + 
	12.0*beta2_12 + z1_2_1 + 2.0*z1_2_2 + 3.0*z1_2_3 + 4.0*z1_2_4 + 5.0*z1_2_5 + 
	beta3_1 + 2.0*beta3_2 + 3.0*beta3_3 + 4.0*beta3_4 + 5.0*beta3_5 + 6.0*beta3_6 + 
	7.0*beta3_7 + 8.0*beta3_8 + 9.0*beta3_9 + z1_3_1 + 2.0*z1_3_2 + 3.0*z1_3_3 + 
	4.0*z1_3_4 + 5.0*z1_3_5 + beta4_1 + 2.0*beta4_2 + 3.0*beta4_3 + 4.0*beta4_4 + 
	5.0*beta4_5 + 6.0*beta4_6 + z1_4_1 + 2.0*z1_4_2 + 3.0*z1_4_3 + 4.0*z1_4_4 + 
	5.0*z1_4_5 + 12.0;
subject to nln172:
	0 >= -(sqrt((m1*n2_1))) - (sqrt((m2*n2_2))) - (sqrt((m3*n2_3))) - 
	(sqrt((m4*n2_4))) + beta1_1 + 2.0*beta1_2 + 3.0*beta1_3 + 4.0*beta1_4 + 
	5.0*beta1_5 + 6.0*beta1_6 + 7.0*beta1_7 + 8.0*beta1_8 + 9.0*beta1_9 + 
	10.0*beta1_10 + 11.0*beta1_11 + 12.0*beta1_12 + 13.0*beta1_13 + 14.0*beta1_14 + 
	15.0*beta1_15 + z2_1_1 + 2.0*z2_1_2 + 3.0*z2_1_3 + 4.0*z2_1_4 + 5.0*z2_1_5 + 
	beta2_1 + 2.0*beta2_2 + 3.0*beta2_3 + 4.0*beta2_4 + 5.0*beta2_5 + 6.0*beta2_6 + 
	7.0*beta2_7 + 8.0*beta2_8 + 9.0*beta2_9 + 10.0*beta2_10 + 11.0*beta2_11 + 
	12.0*beta2_12 + z2_2_1 + 2.0*z2_2_2 + 3.0*z2_2_3 + 4.0*z2_2_4 + 5.0*z2_2_5 + 
	beta3_1 + 2.0*beta3_2 + 3.0*beta3_3 + 4.0*beta3_4 + 5.0*beta3_5 + 6.0*beta3_6 + 
	7.0*beta3_7 + 8.0*beta3_8 + 9.0*beta3_9 + z2_3_1 + 2.0*z2_3_2 + 3.0*z2_3_3 + 
	4.0*z2_3_4 + 5.0*z2_3_5 + beta4_1 + 2.0*beta4_2 + 3.0*beta4_3 + 4.0*beta4_4 + 
	5.0*beta4_5 + 6.0*beta4_6 + z2_4_1 + 2.0*z2_4_2 + 3.0*z2_4_3 + 4.0*z2_4_4 + 
	5.0*z2_4_5 + 6.0;
subject to nln173:
	0 >= -(sqrt((m1*n3_1))) - (sqrt((m2*n3_2))) - (sqrt((m3*n3_3))) - 
	(sqrt((m4*n3_4))) + beta1_1 + 2.0*beta1_2 + 3.0*beta1_3 + 4.0*beta1_4 + 
	5.0*beta1_5 + 6.0*beta1_6 + 7.0*beta1_7 + 8.0*beta1_8 + 9.0*beta1_9 + 
	10.0*beta1_10 + 11.0*beta1_11 + 12.0*beta1_12 + 13.0*beta1_13 + 14.0*beta1_14 + 
	15.0*beta1_15 + z3_1_1 + 2.0*z3_1_2 + 3.0*z3_1_3 + 4.0*z3_1_4 + 5.0*z3_1_5 + 
	beta2_1 + 2.0*beta2_2 + 3.0*beta2_3 + 4.0*beta2_4 + 5.0*beta2_5 + 6.0*beta2_6 + 
	7.0*beta2_7 + 8.0*beta2_8 + 9.0*beta2_9 + 10.0*beta2_10 + 11.0*beta2_11 + 
	12.0*beta2_12 + z3_2_1 + 2.0*z3_2_2 + 3.0*z3_2_3 + 4.0*z3_2_4 + 5.0*z3_2_5 + 
	beta3_1 + 2.0*beta3_2 + 3.0*beta3_3 + 4.0*beta3_4 + 5.0*beta3_5 + 6.0*beta3_6 + 
	7.0*beta3_7 + 8.0*beta3_8 + 9.0*beta3_9 + z3_3_1 + 2.0*z3_3_2 + 3.0*z3_3_3 + 
	4.0*z3_3_4 + 5.0*z3_3_5 + beta4_1 + 2.0*beta4_2 + 3.0*beta4_3 + 4.0*beta4_4 + 
	5.0*beta4_5 + 6.0*beta4_6 + z3_4_1 + 2.0*z3_4_2 + 3.0*z3_4_3 + 4.0*z3_4_4 + 
	5.0*z3_4_5 + 15.0;
subject to nln174:
	0 >= -(sqrt((m1*n4_1))) - (sqrt((m2*n4_2))) - (sqrt((m3*n4_3))) - 
	(sqrt((m4*n4_4))) + beta1_1 + 2.0*beta1_2 + 3.0*beta1_3 + 4.0*beta1_4 + 
	5.0*beta1_5 + 6.0*beta1_6 + 7.0*beta1_7 + 8.0*beta1_8 + 9.0*beta1_9 + 
	10.0*beta1_10 + 11.0*beta1_11 + 12.0*beta1_12 + 13.0*beta1_13 + 14.0*beta1_14 + 
	15.0*beta1_15 + z4_1_1 + 2.0*z4_1_2 + 3.0*z4_1_3 + 4.0*z4_1_4 + beta2_1 + 
	2.0*beta2_2 + 3.0*beta2_3 + 4.0*beta2_4 + 5.0*beta2_5 + 6.0*beta2_6 + 
	7.0*beta2_7 + 8.0*beta2_8 + 9.0*beta2_9 + 10.0*beta2_10 + 11.0*beta2_11 + 
	12.0*beta2_12 + z4_2_1 + 2.0*z4_2_2 + 3.0*z4_2_3 + 4.0*z4_2_4 + beta3_1 + 
	2.0*beta3_2 + 3.0*beta3_3 + 4.0*beta3_4 + 5.0*beta3_5 + 6.0*beta3_6 + 
	7.0*beta3_7 + 8.0*beta3_8 + 9.0*beta3_9 + z4_3_1 + 2.0*z4_3_2 + 3.0*z4_3_3 + 
	4.0*z4_3_4 + beta4_1 + 2.0*beta4_2 + 3.0*beta4_3 + 4.0*beta4_4 + 5.0*beta4_5 + 
	6.0*beta4_6 + z4_4_1 + 2.0*z4_4_2 + 3.0*z4_4_3 + 4.0*z4_4_4 + 9.0;

solve;
	display y1;
	display y2;
	display y3;
	display y4;
	display beta1_1;
	display beta1_2;
	display beta1_3;
	display beta1_4;
	display beta1_5;
	display beta1_6;
	display beta1_7;
	display beta1_8;
	display beta1_9;
	display beta1_10;
	display beta1_11;
	display beta1_12;
	display beta1_13;
	display beta1_14;
	display beta1_15;
	display beta2_1;
	display beta2_2;
	display beta2_3;
	display beta2_4;
	display beta2_5;
	display beta2_6;
	display beta2_7;
	display beta2_8;
	display beta2_9;
	display beta2_10;
	display beta2_11;
	display beta2_12;
	display beta3_1;
	display beta3_2;
	display beta3_3;
	display beta3_4;
	display beta3_5;
	display beta3_6;
	display beta3_7;
	display beta3_8;
	display beta3_9;
	display beta4_1;
	display beta4_2;
	display beta4_3;
	display beta4_4;
	display beta4_5;
	display beta4_6;
	display z1_1_1;
	display z1_1_2;
	display z1_1_3;
	display z1_1_4;
	display z1_1_5;
	display z1_2_1;
	display z1_2_2;
	display z1_2_3;
	display z1_2_4;
	display z1_2_5;
	display z1_3_1;
	display z1_3_2;
	display z1_3_3;
	display z1_3_4;
	display z1_3_5;
	display z1_4_1;
	display z1_4_2;
	display z1_4_3;
	display z1_4_4;
	display z1_4_5;
	display z2_1_1;
	display z2_1_2;
	display z2_1_3;
	display z2_1_4;
	display z2_1_5;
	display z2_2_1;
	display z2_2_2;
	display z2_2_3;
	display z2_2_4;
	display z2_2_5;
	display z2_3_1;
	display z2_3_2;
	display z2_3_3;
	display z2_3_4;
	display z2_3_5;
	display z2_4_1;
	display z2_4_2;
	display z2_4_3;
	display z2_4_4;
	display z2_4_5;
	display z3_1_1;
	display z3_1_2;
	display z3_1_3;
	display z3_1_4;
	display z3_1_5;
	display z3_2_1;
	display z3_2_2;
	display z3_2_3;
	display z3_2_4;
	display z3_2_5;
	display z3_3_1;
	display z3_3_2;
	display z3_3_3;
	display z3_3_4;
	display z3_3_5;
	display z3_4_1;
	display z3_4_2;
	display z3_4_3;
	display z3_4_4;
	display z3_4_5;
	display z4_1_1;
	display z4_1_2;
	display z4_1_3;
	display z4_1_4;
	display z4_2_1;
	display z4_2_2;
	display z4_2_3;
	display z4_2_4;
	display z4_3_1;
	display z4_3_2;
	display z4_3_3;
	display z4_3_4;
	display z4_4_1;
	display z4_4_2;
	display z4_4_3;
	display z4_4_4;
	display m1;
	display m2;
	display m3;
	display m4;
	display n1_1;
	display n1_2;
	display n1_3;
	display n1_4;
	display n2_1;
	display n2_2;
	display n2_3;
	display n2_4;
	display n3_1;
	display n3_2;
	display n3_3;
	display n3_4;
	display n4_1;
	display n4_2;
	display n4_3;
	display n4_4;
display obj;
