#
#**************************
# SET UP THE INITIAL DATA *
#**************************
#   Problem :
#   *********
#
#   A nonlinear formulation of a 0-1 combinatorial problem
#   which can be described as grouping and clustering.
#   This version of the problem consists in grouping 20 cans
#   into 5 groups and the objective is to obtain groups that
#   are as homogeneous as possible with respect to the characteristics
#   of the cans.
#   Note that the initial point is feasible.
#   This problem is the contribution of a LANCELOT user in order to fulfill
#   the package licence agreement.
#   Source :
#   D. Tuyttens, M. Pirlot, J. Teghem, E. Trauwaert and B. Liegeois,
#   "Homogeneous grouping of nuclear fuel cans through simulated
#   annealing and tabu search",
#   Annals of Operations Research 50 (1994) 575-607.
#   SIF input: D. Tuyttens (FPMs, Mons, Belgium) Dec 1994.
#   classification SQR2-MN-100-125
#   Problem parameters : Number of cans, groups and characteristics
#   The cans and their characteristics
#   First characteristic
#   Second characteristic
#   Average characteristics
#   Variables of the problem
#   X(j,k) = 1  can j belongs to group k
#          = 0  otherwise
#   Objective function and Constraints
#   Linear constraints
#   Non-Linear constraints
# The constants of the problem.
# The bounds on the variables
# The starting point
#  A feasible solution
#  This solution consists into setting the first 4 cans in 
#  group 1 and  the next 4 ones in the next group and so on.
# The nonlinear elements functions
# The group functions
# Solution
	param lots := 5;
	param pots := 20;
	param car := 2;
	param t1_1 := 4.37;
	param t2_1 := 4.56;
	param t3_1 := 4.26;
	param t4_1 := 4.56;
	param t5_1 := 4.3;
	param t6_1 := 4.46;
	param t7_1 := 3.84;
	param t8_1 := 4.57;
	param t9_1 := 4.26;
	param t10_1 := 4.37;
	param t11_1 := 3.49;
	param t12_1 := 4.43;
	param t13_1 := 4.48;
	param t14_1 := 4.01;
	param t15_1 := 4.29;
	param t16_1 := 4.42;
	param t17_1 := 4.23;
	param t18_1 := 4.42;
	param t19_1 := 4.23;
	param t20_1 := 3.49;
	param t1_2 := 5.23;
	param t2_2 := 5.74;
	param t3_2 := 4.93;
	param t4_2 := 5.74;
	param t5_2 := 5.19;
	param t6_2 := 5.46;
	param t7_2 := 4.65;
	param t8_2 := 5.27;
	param t9_2 := 5.57;
	param t10_2 := 5.12;
	param t11_2 := 5.73;
	param t12_2 := 5.45;
	param t13_2 := 5.42;
	param t14_2 := 4.05;
	param t15_2 := 4.26;
	param t16_2 := 4.58;
	param t17_2 := 3.94;
	param t18_2 := 4.18;
	param t19_2 := 4.18;
	param t20_2 := 5.89;
	param xm1 := 17.1;
	param xm2 := 20.1;
	param potlot := 4.0;

	var x1_1 >= 0.0 ,  <= 1.0 ,  := 1.0;
	var x1_2 >= 0.0 ,  <= 1.0;
	var x1_3 >= 0.0 ,  <= 1.0;
	var x1_4 >= 0.0 ,  <= 1.0;
	var x1_5 >= 0.0 ,  <= 1.0;
	var x2_1 >= 0.0 ,  <= 1.0 ,  := 1.0;
	var x2_2 >= 0.0 ,  <= 1.0;
	var x2_3 >= 0.0 ,  <= 1.0;
	var x2_4 >= 0.0 ,  <= 1.0;
	var x2_5 >= 0.0 ,  <= 1.0;
	var x3_1 >= 0.0 ,  <= 1.0 ,  := 1.0;
	var x3_2 >= 0.0 ,  <= 1.0;
	var x3_3 >= 0.0 ,  <= 1.0;
	var x3_4 >= 0.0 ,  <= 1.0;
	var x3_5 >= 0.0 ,  <= 1.0;
	var x4_1 >= 0.0 ,  <= 1.0 ,  := 1.0;
	var x4_2 >= 0.0 ,  <= 1.0;
	var x4_3 >= 0.0 ,  <= 1.0;
	var x4_4 >= 0.0 ,  <= 1.0;
	var x4_5 >= 0.0 ,  <= 1.0;
	var x5_1 >= 0.0 ,  <= 1.0;
	var x5_2 >= 0.0 ,  <= 1.0 ,  := 1.0;
	var x5_3 >= 0.0 ,  <= 1.0;
	var x5_4 >= 0.0 ,  <= 1.0;
	var x5_5 >= 0.0 ,  <= 1.0;
	var x6_1 >= 0.0 ,  <= 1.0;
	var x6_2 >= 0.0 ,  <= 1.0 ,  := 1.0;
	var x6_3 >= 0.0 ,  <= 1.0;
	var x6_4 >= 0.0 ,  <= 1.0;
	var x6_5 >= 0.0 ,  <= 1.0;
	var x7_1 >= 0.0 ,  <= 1.0;
	var x7_2 >= 0.0 ,  <= 1.0 ,  := 1.0;
	var x7_3 >= 0.0 ,  <= 1.0;
	var x7_4 >= 0.0 ,  <= 1.0;
	var x7_5 >= 0.0 ,  <= 1.0;
	var x8_1 >= 0.0 ,  <= 1.0;
	var x8_2 >= 0.0 ,  <= 1.0 ,  := 1.0;
	var x8_3 >= 0.0 ,  <= 1.0;
	var x8_4 >= 0.0 ,  <= 1.0;
	var x8_5 >= 0.0 ,  <= 1.0;
	var x9_1 >= 0.0 ,  <= 1.0;
	var x9_2 >= 0.0 ,  <= 1.0;
	var x9_3 >= 0.0 ,  <= 1.0 ,  := 1.0;
	var x9_4 >= 0.0 ,  <= 1.0;
	var x9_5 >= 0.0 ,  <= 1.0;
	var x10_1 >= 0.0 ,  <= 1.0;
	var x10_2 >= 0.0 ,  <= 1.0;
	var x10_3 >= 0.0 ,  <= 1.0 ,  := 1.0;
	var x10_4 >= 0.0 ,  <= 1.0;
	var x10_5 >= 0.0 ,  <= 1.0;
	var x11_1 >= 0.0 ,  <= 1.0;
	var x11_2 >= 0.0 ,  <= 1.0;
	var x11_3 >= 0.0 ,  <= 1.0 ,  := 1.0;
	var x11_4 >= 0.0 ,  <= 1.0;
	var x11_5 >= 0.0 ,  <= 1.0;
	var x12_1 >= 0.0 ,  <= 1.0;
	var x12_2 >= 0.0 ,  <= 1.0;
	var x12_3 >= 0.0 ,  <= 1.0 ,  := 1.0;
	var x12_4 >= 0.0 ,  <= 1.0;
	var x12_5 >= 0.0 ,  <= 1.0;
	var x13_1 >= 0.0 ,  <= 1.0;
	var x13_2 >= 0.0 ,  <= 1.0;
	var x13_3 >= 0.0 ,  <= 1.0;
	var x13_4 >= 0.0 ,  <= 1.0 ,  := 1.0;
	var x13_5 >= 0.0 ,  <= 1.0;
	var x14_1 >= 0.0 ,  <= 1.0;
	var x14_2 >= 0.0 ,  <= 1.0;
	var x14_3 >= 0.0 ,  <= 1.0;
	var x14_4 >= 0.0 ,  <= 1.0 ,  := 1.0;
	var x14_5 >= 0.0 ,  <= 1.0;
	var x15_1 >= 0.0 ,  <= 1.0;
	var x15_2 >= 0.0 ,  <= 1.0;
	var x15_3 >= 0.0 ,  <= 1.0;
	var x15_4 >= 0.0 ,  <= 1.0 ,  := 1.0;
	var x15_5 >= 0.0 ,  <= 1.0;
	var x16_1 >= 0.0 ,  <= 1.0;
	var x16_2 >= 0.0 ,  <= 1.0;
	var x16_3 >= 0.0 ,  <= 1.0;
	var x16_4 >= 0.0 ,  <= 1.0 ,  := 1.0;
	var x16_5 >= 0.0 ,  <= 1.0;
	var x17_1 >= 0.0 ,  <= 1.0;
	var x17_2 >= 0.0 ,  <= 1.0;
	var x17_3 >= 0.0 ,  <= 1.0;
	var x17_4 >= 0.0 ,  <= 1.0;
	var x17_5 >= 0.0 ,  <= 1.0 ,  := 1.0;
	var x18_1 >= 0.0 ,  <= 1.0;
	var x18_2 >= 0.0 ,  <= 1.0;
	var x18_3 >= 0.0 ,  <= 1.0;
	var x18_4 >= 0.0 ,  <= 1.0;
	var x18_5 >= 0.0 ,  <= 1.0 ,  := 1.0;
	var x19_1 >= 0.0 ,  <= 1.0;
	var x19_2 >= 0.0 ,  <= 1.0;
	var x19_3 >= 0.0 ,  <= 1.0;
	var x19_4 >= 0.0 ,  <= 1.0;
	var x19_5 >= 0.0 ,  <= 1.0 ,  := 1.0;
	var x20_1 >= 0.0 ,  <= 1.0;
	var x20_2 >= 0.0 ,  <= 1.0;
	var x20_3 >= 0.0 ,  <= 1.0;
	var x20_4 >= 0.0 ,  <= 1.0;
	var x20_5 >= 0.0 ,  <= 1.0 ,  := 1.0;

minimize obj:
	(4.37*x1_1 + 4.56*x2_1 + 4.26*x3_1 + 4.56*x4_1 + 4.3*x5_1 + 4.46*x6_1 + 
	3.84*x7_1 + 4.57*x8_1 + 4.26*x9_1 + 4.37*x10_1 + 3.49*x11_1 + 4.43*x12_1 + 
	4.48*x13_1 + 4.01*x14_1 + 4.29*x15_1 + 4.42*x16_1 + 4.23*x17_1 + 4.42*x18_1 + 
	4.23*x19_1 + 3.49*x20_1 - 17.1)*(4.37*x1_1 + 4.56*x2_1 + 4.26*x3_1 + 4.56*x4_1 
	+ 4.3*x5_1 + 4.46*x6_1 + 3.84*x7_1 + 4.57*x8_1 + 4.26*x9_1 + 4.37*x10_1 + 
	3.49*x11_1 + 4.43*x12_1 + 4.48*x13_1 + 4.01*x14_1 + 4.29*x15_1 + 4.42*x16_1 + 
	4.23*x17_1 + 4.42*x18_1 + 4.23*x19_1 + 3.49*x20_1 - 17.1) + (5.23*x1_1 + 
	5.74*x2_1 + 4.93*x3_1 + 5.74*x4_1 + 5.19*x5_1 + 5.46*x6_1 + 4.65*x7_1 + 
	5.27*x8_1 + 5.57*x9_1 + 5.12*x10_1 + 5.73*x11_1 + 5.45*x12_1 + 5.42*x13_1 + 
	4.05*x14_1 + 4.26*x15_1 + 4.58*x16_1 + 3.94*x17_1 + 4.18*x18_1 + 4.18*x19_1 + 
	5.89*x20_1 - 20.1)*(5.23*x1_1 + 5.74*x2_1 + 4.93*x3_1 + 5.74*x4_1 + 5.19*x5_1 + 
	5.46*x6_1 + 4.65*x7_1 + 5.27*x8_1 + 5.57*x9_1 + 5.12*x10_1 + 5.73*x11_1 + 
	5.45*x12_1 + 5.42*x13_1 + 4.05*x14_1 + 4.26*x15_1 + 4.58*x16_1 + 3.94*x17_1 + 
	4.18*x18_1 + 4.18*x19_1 + 5.89*x20_1 - 20.1) + (4.37*x1_2 + 4.56*x2_2 + 
	4.26*x3_2 + 4.56*x4_2 + 4.3*x5_2 + 4.46*x6_2 + 3.84*x7_2 + 4.57*x8_2 + 
	4.26*x9_2 + 4.37*x10_2 + 3.49*x11_2 + 4.43*x12_2 + 4.48*x13_2 + 4.01*x14_2 + 
	4.29*x15_2 + 4.42*x16_2 + 4.23*x17_2 + 4.42*x18_2 + 4.23*x19_2 + 3.49*x20_2 - 
	17.1)*(4.37*x1_2 + 4.56*x2_2 + 4.26*x3_2 + 4.56*x4_2 + 4.3*x5_2 + 4.46*x6_2 + 
	3.84*x7_2 + 4.57*x8_2 + 4.26*x9_2 + 4.37*x10_2 + 3.49*x11_2 + 4.43*x12_2 + 
	4.48*x13_2 + 4.01*x14_2 + 4.29*x15_2 + 4.42*x16_2 + 4.23*x17_2 + 4.42*x18_2 + 
	4.23*x19_2 + 3.49*x20_2 - 17.1) + (5.23*x1_2 + 5.74*x2_2 + 4.93*x3_2 + 
	5.74*x4_2 + 5.19*x5_2 + 5.46*x6_2 + 4.65*x7_2 + 5.27*x8_2 + 5.57*x9_2 + 
	5.12*x10_2 + 5.73*x11_2 + 5.45*x12_2 + 5.42*x13_2 + 4.05*x14_2 + 4.26*x15_2 + 
	4.58*x16_2 + 3.94*x17_2 + 4.18*x18_2 + 4.18*x19_2 + 5.89*x20_2 - 
	20.1)*(5.23*x1_2 + 5.74*x2_2 + 4.93*x3_2 + 5.74*x4_2 + 5.19*x5_2 + 5.46*x6_2 + 
	4.65*x7_2 + 5.27*x8_2 + 5.57*x9_2 + 5.12*x10_2 + 5.73*x11_2 + 5.45*x12_2 + 
	5.42*x13_2 + 4.05*x14_2 + 4.26*x15_2 + 4.58*x16_2 + 3.94*x17_2 + 4.18*x18_2 + 
	4.18*x19_2 + 5.89*x20_2 - 20.1) + (4.37*x1_3 + 4.56*x2_3 + 4.26*x3_3 + 
	4.56*x4_3 + 4.3*x5_3 + 4.46*x6_3 + 3.84*x7_3 + 4.57*x8_3 + 4.26*x9_3 + 
	4.37*x10_3 + 3.49*x11_3 + 4.43*x12_3 + 4.48*x13_3 + 4.01*x14_3 + 4.29*x15_3 + 
	4.42*x16_3 + 4.23*x17_3 + 4.42*x18_3 + 4.23*x19_3 + 3.49*x20_3 - 
	17.1)*(4.37*x1_3 + 4.56*x2_3 + 4.26*x3_3 + 4.56*x4_3 + 4.3*x5_3 + 4.46*x6_3 + 
	3.84*x7_3 + 4.57*x8_3 + 4.26*x9_3 + 4.37*x10_3 + 3.49*x11_3 + 4.43*x12_3 + 
	4.48*x13_3 + 4.01*x14_3 + 4.29*x15_3 + 4.42*x16_3 + 4.23*x17_3 + 4.42*x18_3 + 
	4.23*x19_3 + 3.49*x20_3 - 17.1) + (5.23*x1_3 + 5.74*x2_3 + 4.93*x3_3 + 
	5.74*x4_3 + 5.19*x5_3 + 5.46*x6_3 + 4.65*x7_3 + 5.27*x8_3 + 5.57*x9_3 + 
	5.12*x10_3 + 5.73*x11_3 + 5.45*x12_3 + 5.42*x13_3 + 4.05*x14_3 + 4.26*x15_3 + 
	4.58*x16_3 + 3.94*x17_3 + 4.18*x18_3 + 4.18*x19_3 + 5.89*x20_3 - 
	20.1)*(5.23*x1_3 + 5.74*x2_3 + 4.93*x3_3 + 5.74*x4_3 + 5.19*x5_3 + 5.46*x6_3 + 
	4.65*x7_3 + 5.27*x8_3 + 5.57*x9_3 + 5.12*x10_3 + 5.73*x11_3 + 5.45*x12_3 + 
	5.42*x13_3 + 4.05*x14_3 + 4.26*x15_3 + 4.58*x16_3 + 3.94*x17_3 + 4.18*x18_3 + 
	4.18*x19_3 + 5.89*x20_3 - 20.1) + (4.37*x1_4 + 4.56*x2_4 + 4.26*x3_4 + 
	4.56*x4_4 + 4.3*x5_4 + 4.46*x6_4 + 3.84*x7_4 + 4.57*x8_4 + 4.26*x9_4 + 
	4.37*x10_4 + 3.49*x11_4 + 4.43*x12_4 + 4.48*x13_4 + 4.01*x14_4 + 4.29*x15_4 + 
	4.42*x16_4 + 4.23*x17_4 + 4.42*x18_4 + 4.23*x19_4 + 3.49*x20_4 - 
	17.1)*(4.37*x1_4 + 4.56*x2_4 + 4.26*x3_4 + 4.56*x4_4 + 4.3*x5_4 + 4.46*x6_4 + 
	3.84*x7_4 + 4.57*x8_4 + 4.26*x9_4 + 4.37*x10_4 + 3.49*x11_4 + 4.43*x12_4 + 
	4.48*x13_4 + 4.01*x14_4 + 4.29*x15_4 + 4.42*x16_4 + 4.23*x17_4 + 4.42*x18_4 + 
	4.23*x19_4 + 3.49*x20_4 - 17.1) + (5.23*x1_4 + 5.74*x2_4 + 4.93*x3_4 + 
	5.74*x4_4 + 5.19*x5_4 + 5.46*x6_4 + 4.65*x7_4 + 5.27*x8_4 + 5.57*x9_4 + 
	5.12*x10_4 + 5.73*x11_4 + 5.45*x12_4 + 5.42*x13_4 + 4.05*x14_4 + 4.26*x15_4 + 
	4.58*x16_4 + 3.94*x17_4 + 4.18*x18_4 + 4.18*x19_4 + 5.89*x20_4 - 
	20.1)*(5.23*x1_4 + 5.74*x2_4 + 4.93*x3_4 + 5.74*x4_4 + 5.19*x5_4 + 5.46*x6_4 + 
	4.65*x7_4 + 5.27*x8_4 + 5.57*x9_4 + 5.12*x10_4 + 5.73*x11_4 + 5.45*x12_4 + 
	5.42*x13_4 + 4.05*x14_4 + 4.26*x15_4 + 4.58*x16_4 + 3.94*x17_4 + 4.18*x18_4 + 
	4.18*x19_4 + 5.89*x20_4 - 20.1) + (4.37*x1_5 + 4.56*x2_5 + 4.26*x3_5 + 
	4.56*x4_5 + 4.3*x5_5 + 4.46*x6_5 + 3.84*x7_5 + 4.57*x8_5 + 4.26*x9_5 + 
	4.37*x10_5 + 3.49*x11_5 + 4.43*x12_5 + 4.48*x13_5 + 4.01*x14_5 + 4.29*x15_5 + 
	4.42*x16_5 + 4.23*x17_5 + 4.42*x18_5 + 4.23*x19_5 + 3.49*x20_5 - 
	17.1)*(4.37*x1_5 + 4.56*x2_5 + 4.26*x3_5 + 4.56*x4_5 + 4.3*x5_5 + 4.46*x6_5 + 
	3.84*x7_5 + 4.57*x8_5 + 4.26*x9_5 + 4.37*x10_5 + 3.49*x11_5 + 4.43*x12_5 + 
	4.48*x13_5 + 4.01*x14_5 + 4.29*x15_5 + 4.42*x16_5 + 4.23*x17_5 + 4.42*x18_5 + 
	4.23*x19_5 + 3.49*x20_5 - 17.1) + (5.23*x1_5 + 5.74*x2_5 + 4.93*x3_5 + 
	5.74*x4_5 + 5.19*x5_5 + 5.46*x6_5 + 4.65*x7_5 + 5.27*x8_5 + 5.57*x9_5 + 
	5.12*x10_5 + 5.73*x11_5 + 5.45*x12_5 + 5.42*x13_5 + 4.05*x14_5 + 4.26*x15_5 + 
	4.58*x16_5 + 3.94*x17_5 + 4.18*x18_5 + 4.18*x19_5 + 5.89*x20_5 - 
	20.1)*(5.23*x1_5 + 5.74*x2_5 + 4.93*x3_5 + 5.74*x4_5 + 5.19*x5_5 + 5.46*x6_5 + 
	4.65*x7_5 + 5.27*x8_5 + 5.57*x9_5 + 5.12*x10_5 + 5.73*x11_5 + 5.45*x12_5 + 
	5.42*x13_5 + 4.05*x14_5 + 4.26*x15_5 + 4.58*x16_5 + 3.94*x17_5 + 4.18*x18_5 + 
	4.18*x19_5 + 5.89*x20_5 - 20.1);

subject to cl11:
	x1_1 + x2_1 + x3_1 + x4_1 + x5_1 + x6_1 + x7_1 + x8_1 + x9_1 + x10_1 + x11_1 + 
	x12_1 + x13_1 + x14_1 + x15_1 + x16_1 + x17_1 + x18_1 + x19_1 + x20_1 - 4.0 = 0;
subject to cl12:
	x1_2 + x2_2 + x3_2 + x4_2 + x5_2 + x6_2 + x7_2 + x8_2 + x9_2 + x10_2 + x11_2 + 
	x12_2 + x13_2 + x14_2 + x15_2 + x16_2 + x17_2 + x18_2 + x19_2 + x20_2 - 4.0 = 0;
subject to cl13:
	x1_3 + x2_3 + x3_3 + x4_3 + x5_3 + x6_3 + x7_3 + x8_3 + x9_3 + x10_3 + x11_3 + 
	x12_3 + x13_3 + x14_3 + x15_3 + x16_3 + x17_3 + x18_3 + x19_3 + x20_3 - 4.0 = 0;
subject to cl14:
	x1_4 + x2_4 + x3_4 + x4_4 + x5_4 + x6_4 + x7_4 + x8_4 + x9_4 + x10_4 + x11_4 + 
	x12_4 + x13_4 + x14_4 + x15_4 + x16_4 + x17_4 + x18_4 + x19_4 + x20_4 - 4.0 = 0;
subject to cl15:
	x1_5 + x2_5 + x3_5 + x4_5 + x5_5 + x6_5 + x7_5 + x8_5 + x9_5 + x10_5 + x11_5 + 
	x12_5 + x13_5 + x14_5 + x15_5 + x16_5 + x17_5 + x18_5 + x19_5 + x20_5 - 4.0 = 0;
subject to cl21:
	x1_1 + x1_2 + x1_3 + x1_4 + x1_5 - 1.0 = 0;
subject to cl22:
	x2_1 + x2_2 + x2_3 + x2_4 + x2_5 - 1.0 = 0;
subject to cl23:
	x3_1 + x3_2 + x3_3 + x3_4 + x3_5 - 1.0 = 0;
subject to cl24:
	x4_1 + x4_2 + x4_3 + x4_4 + x4_5 - 1.0 = 0;
subject to cl25:
	x5_1 + x5_2 + x5_3 + x5_4 + x5_5 - 1.0 = 0;
subject to cl26:
	x6_1 + x6_2 + x6_3 + x6_4 + x6_5 - 1.0 = 0;
subject to cl27:
	x7_1 + x7_2 + x7_3 + x7_4 + x7_5 - 1.0 = 0;
subject to cl28:
	x8_1 + x8_2 + x8_3 + x8_4 + x8_5 - 1.0 = 0;
subject to cl29:
	x9_1 + x9_2 + x9_3 + x9_4 + x9_5 - 1.0 = 0;
subject to cl210:
	x10_1 + x10_2 + x10_3 + x10_4 + x10_5 - 1.0 = 0;
subject to cl211:
	x11_1 + x11_2 + x11_3 + x11_4 + x11_5 - 1.0 = 0;
subject to cl212:
	x12_1 + x12_2 + x12_3 + x12_4 + x12_5 - 1.0 = 0;
subject to cl213:
	x13_1 + x13_2 + x13_3 + x13_4 + x13_5 - 1.0 = 0;
subject to cl214:
	x14_1 + x14_2 + x14_3 + x14_4 + x14_5 - 1.0 = 0;
subject to cl215:
	x15_1 + x15_2 + x15_3 + x15_4 + x15_5 - 1.0 = 0;
subject to cl216:
	x16_1 + x16_2 + x16_3 + x16_4 + x16_5 - 1.0 = 0;
subject to cl217:
	x17_1 + x17_2 + x17_3 + x17_4 + x17_5 - 1.0 = 0;
subject to cl218:
	x18_1 + x18_2 + x18_3 + x18_4 + x18_5 - 1.0 = 0;
subject to cl219:
	x19_1 + x19_2 + x19_3 + x19_4 + x19_5 - 1.0 = 0;
subject to cl220:
	x20_1 + x20_2 + x20_3 + x20_4 + x20_5 - 1.0 = 0;
subject to cnl1_1:
	x1_1 * (x1_1 - 1.0 )  = 0;
subject to cnl1_2:
	x1_2 * (x1_2 - 1.0 )  = 0;
subject to cnl1_3:
	x1_3 * (x1_3 - 1.0 )  = 0;
subject to cnl1_4:
	x1_4 * (x1_4 - 1.0 )  = 0;
subject to cnl1_5:
	x1_5 * (x1_5 - 1.0 )  = 0;
subject to cnl2_1:
	x2_1 * (x2_1 - 1.0 )  = 0;
subject to cnl2_2:
	x2_2 * (x2_2 - 1.0 )  = 0;
subject to cnl2_3:
	x2_3 * (x2_3 - 1.0 )  = 0;
subject to cnl2_4:
	x2_4 * (x2_4 - 1.0 )  = 0;
subject to cnl2_5:
	x2_5 * (x2_5 - 1.0 )  = 0;
subject to cnl3_1:
	x3_1 * (x3_1 - 1.0 )  = 0;
subject to cnl3_2:
	x3_2 * (x3_2 - 1.0 )  = 0;
subject to cnl3_3:
	x3_3 * (x3_3 - 1.0 )  = 0;
subject to cnl3_4:
	x3_4 * (x3_4 - 1.0 )  = 0;
subject to cnl3_5:
	x3_5 * (x3_5 - 1.0 )  = 0;
subject to cnl4_1:
	x4_1 * (x4_1 - 1.0 )  = 0;
subject to cnl4_2:
	x4_2 * (x4_2 - 1.0 )  = 0;
subject to cnl4_3:
	x4_3 * (x4_3 - 1.0 )  = 0;
subject to cnl4_4:
	x4_4 * (x4_4 - 1.0 )  = 0;
subject to cnl4_5:
	x4_5 * (x4_5 - 1.0 )  = 0;
subject to cnl5_1:
	x5_1 * (x5_1 - 1.0 )  = 0;
subject to cnl5_2:
	x5_2 * (x5_2 - 1.0 )  = 0;
subject to cnl5_3:
	x5_3 * (x5_3 - 1.0 )  = 0;
subject to cnl5_4:
	x5_4 * (x5_4 - 1.0 )  = 0;
subject to cnl5_5:
	x5_5 * (x5_5 - 1.0 )  = 0;
subject to cnl6_1:
	x6_1 * (x6_1 - 1.0 )  = 0;
subject to cnl6_2:
	x6_2 * (x6_2 - 1.0 )  = 0;
subject to cnl6_3:
	x6_3 * (x6_3 - 1.0 )  = 0;
subject to cnl6_4:
	x6_4 * (x6_4 - 1.0 )  = 0;
subject to cnl6_5:
	x6_5 * (x6_5 - 1.0 )  = 0;
subject to cnl7_1:
	x7_1 * (x7_1 - 1.0 )  = 0;
subject to cnl7_2:
	x7_2 * (x7_2 - 1.0 )  = 0;
subject to cnl7_3:
	x7_3 * (x7_3 - 1.0 )  = 0;
subject to cnl7_4:
	x7_4 * (x7_4 - 1.0 )  = 0;
subject to cnl7_5:
	x7_5 * (x7_5 - 1.0 )  = 0;
subject to cnl8_1:
	x8_1 * (x8_1 - 1.0 )  = 0;
subject to cnl8_2:
	x8_2 * (x8_2 - 1.0 )  = 0;
subject to cnl8_3:
	x8_3 * (x8_3 - 1.0 )  = 0;
subject to cnl8_4:
	x8_4 * (x8_4 - 1.0 )  = 0;
subject to cnl8_5:
	x8_5 * (x8_5 - 1.0 )  = 0;
subject to cnl9_1:
	x9_1 * (x9_1 - 1.0 )  = 0;
subject to cnl9_2:
	x9_2 * (x9_2 - 1.0 )  = 0;
subject to cnl9_3:
	x9_3 * (x9_3 - 1.0 )  = 0;
subject to cnl9_4:
	x9_4 * (x9_4 - 1.0 )  = 0;
subject to cnl9_5:
	x9_5 * (x9_5 - 1.0 )  = 0;
subject to cnl10_1:
	x10_1 * (x10_1 - 1.0 )  = 0;
subject to cnl10_2:
	x10_2 * (x10_2 - 1.0 )  = 0;
subject to cnl10_3:
	x10_3 * (x10_3 - 1.0 )  = 0;
subject to cnl10_4:
	x10_4 * (x10_4 - 1.0 )  = 0;
subject to cnl10_5:
	x10_5 * (x10_5 - 1.0 )  = 0;
subject to cnl11_1:
	x11_1 * (x11_1 - 1.0 )  = 0;
subject to cnl11_2:
	x11_2 * (x11_2 - 1.0 )  = 0;
subject to cnl11_3:
	x11_3 * (x11_3 - 1.0 )  = 0;
subject to cnl11_4:
	x11_4 * (x11_4 - 1.0 )  = 0;
subject to cnl11_5:
	x11_5 * (x11_5 - 1.0 )  = 0;
subject to cnl12_1:
	x12_1 * (x12_1 - 1.0 )  = 0;
subject to cnl12_2:
	x12_2 * (x12_2 - 1.0 )  = 0;
subject to cnl12_3:
	x12_3 * (x12_3 - 1.0 )  = 0;
subject to cnl12_4:
	x12_4 * (x12_4 - 1.0 )  = 0;
subject to cnl12_5:
	x12_5 * (x12_5 - 1.0 )  = 0;
subject to cnl13_1:
	x13_1 * (x13_1 - 1.0 )  = 0;
subject to cnl13_2:
	x13_2 * (x13_2 - 1.0 )  = 0;
subject to cnl13_3:
	x13_3 * (x13_3 - 1.0 )  = 0;
subject to cnl13_4:
	x13_4 * (x13_4 - 1.0 )  = 0;
subject to cnl13_5:
	x13_5 * (x13_5 - 1.0 )  = 0;
subject to cnl14_1:
	x14_1 * (x14_1 - 1.0 )  = 0;
subject to cnl14_2:
	x14_2 * (x14_2 - 1.0 )  = 0;
subject to cnl14_3:
	x14_3 * (x14_3 - 1.0 )  = 0;
subject to cnl14_4:
	x14_4 * (x14_4 - 1.0 )  = 0;
subject to cnl14_5:
	x14_5 * (x14_5 - 1.0 )  = 0;
subject to cnl15_1:
	x15_1 * (x15_1 - 1.0 )  = 0;
subject to cnl15_2:
	x15_2 * (x15_2 - 1.0 )  = 0;
subject to cnl15_3:
	x15_3 * (x15_3 - 1.0 )  = 0;
subject to cnl15_4:
	x15_4 * (x15_4 - 1.0 )  = 0;
subject to cnl15_5:
	x15_5 * (x15_5 - 1.0 )  = 0;
subject to cnl16_1:
	x16_1 * (x16_1 - 1.0 )  = 0;
subject to cnl16_2:
	x16_2 * (x16_2 - 1.0 )  = 0;
subject to cnl16_3:
	x16_3 * (x16_3 - 1.0 )  = 0;
subject to cnl16_4:
	x16_4 * (x16_4 - 1.0 )  = 0;
subject to cnl16_5:
	x16_5 * (x16_5 - 1.0 )  = 0;
subject to cnl17_1:
	x17_1 * (x17_1 - 1.0 )  = 0;
subject to cnl17_2:
	x17_2 * (x17_2 - 1.0 )  = 0;
subject to cnl17_3:
	x17_3 * (x17_3 - 1.0 )  = 0;
subject to cnl17_4:
	x17_4 * (x17_4 - 1.0 )  = 0;
subject to cnl17_5:
	x17_5 * (x17_5 - 1.0 )  = 0;
subject to cnl18_1:
	x18_1 * (x18_1 - 1.0 )  = 0;
subject to cnl18_2:
	x18_2 * (x18_2 - 1.0 )  = 0;
subject to cnl18_3:
	x18_3 * (x18_3 - 1.0 )  = 0;
subject to cnl18_4:
	x18_4 * (x18_4 - 1.0 )  = 0;
subject to cnl18_5:
	x18_5 * (x18_5 - 1.0 )  = 0;
subject to cnl19_1:
	x19_1 * (x19_1 - 1.0 )  = 0;
subject to cnl19_2:
	x19_2 * (x19_2 - 1.0 )  = 0;
subject to cnl19_3:
	x19_3 * (x19_3 - 1.0 )  = 0;
subject to cnl19_4:
	x19_4 * (x19_4 - 1.0 )  = 0;
subject to cnl19_5:
	x19_5 * (x19_5 - 1.0 )  = 0;
subject to cnl20_1:
	x20_1 * (x20_1 - 1.0 )  = 0;
subject to cnl20_2:
	x20_2 * (x20_2 - 1.0 )  = 0;
subject to cnl20_3:
	x20_3 * (x20_3 - 1.0 )  = 0;
subject to cnl20_4:
	x20_4 * (x20_4 - 1.0 )  = 0;
subject to cnl20_5:
	x20_5 * (x20_5 - 1.0 )  = 0;

solve;
	display x1_1;
	display x1_2;
	display x1_3;
	display x1_4;
	display x1_5;
	display x2_1;
	display x2_2;
	display x2_3;
	display x2_4;
	display x2_5;
	display x3_1;
	display x3_2;
	display x3_3;
	display x3_4;
	display x3_5;
	display x4_1;
	display x4_2;
	display x4_3;
	display x4_4;
	display x4_5;
	display x5_1;
	display x5_2;
	display x5_3;
	display x5_4;
	display x5_5;
	display x6_1;
	display x6_2;
	display x6_3;
	display x6_4;
	display x6_5;
	display x7_1;
	display x7_2;
	display x7_3;
	display x7_4;
	display x7_5;
	display x8_1;
	display x8_2;
	display x8_3;
	display x8_4;
	display x8_5;
	display x9_1;
	display x9_2;
	display x9_3;
	display x9_4;
	display x9_5;
	display x10_1;
	display x10_2;
	display x10_3;
	display x10_4;
	display x10_5;
	display x11_1;
	display x11_2;
	display x11_3;
	display x11_4;
	display x11_5;
	display x12_1;
	display x12_2;
	display x12_3;
	display x12_4;
	display x12_5;
	display x13_1;
	display x13_2;
	display x13_3;
	display x13_4;
	display x13_5;
	display x14_1;
	display x14_2;
	display x14_3;
	display x14_4;
	display x14_5;
	display x15_1;
	display x15_2;
	display x15_3;
	display x15_4;
	display x15_5;
	display x16_1;
	display x16_2;
	display x16_3;
	display x16_4;
	display x16_5;
	display x17_1;
	display x17_2;
	display x17_3;
	display x17_4;
	display x17_5;
	display x18_1;
	display x18_2;
	display x18_3;
	display x18_4;
	display x18_5;
	display x19_1;
	display x19_2;
	display x19_3;
	display x19_4;
	display x19_5;
	display x20_1;
	display x20_2;
	display x20_3;
	display x20_4;
	display x20_5;
display obj;
