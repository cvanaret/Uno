#
#**************************
# SET UP THE INITIAL DATA *
#**************************
#   Problem :
#   *********
#   An expanded system formulation of a 3-D PDE system.
#   A nine-point discretization of Laplace's equation in a 
#   rectangular domain may be expressed in the form 
#         - M v = b, 
#   where M = sum a_i a_i^T. Letting A = (a_1 .... a_m), 
#   this system may be expanded as
#          ( I   A^T ) (x) = (0),
#          ( A    0  ) (v)   (b)
#   which is then equivalentto solving the EQP 
#   minimize 1/2 || x ||_2^2   s.t.    A x = b
#   In this variant, we replace the leading I block in the
#   above formulation with a zero-one diagonal matrix D.
#   This corresponds to certain boundary conditions.
#   The resulting QP is thus convex but not strictly convex.
#   SIF input: Nick Gould, February 1994
#   classification QLR2-AN-V-V
#   Number of nodes in x direction
#IE NX                  3
#   Number of nodes in y direction
#IE NY                  3
#   Number of nodes in z direction
#IE NZ                  3
#   Other useful parameters
#  It is easier to describe this problem by columns.
# objective function
# constraints
# objective function terms
# constraints : central constraints
	param nx := 10;
	param ny := 10;
	param nz := 10;
	param xp := 1 + (10);
	param xm := -1 + (10);
	param yp := 1 + (10);
	param ym := -1 + (10);
	param zp := 1 + (10);
	param zm := -1 + (10);
	param m := 10;
	param n := 10;
	param p := 10;
	param kp := 1 + (9);
	param jp := 1 + (9);
	param ip := 1 + (9);

	var x1_1_1;
	var y1_1_1;
	var z1_1_1;
	var x2_1_1;
	var y2_1_1;
	var z2_1_1;
	var x3_1_1;
	var y3_1_1;
	var z3_1_1;
	var x4_1_1;
	var y4_1_1;
	var z4_1_1;
	var x5_1_1;
	var y5_1_1;
	var z5_1_1;
	var x6_1_1;
	var y6_1_1;
	var z6_1_1;
	var x7_1_1;
	var y7_1_1;
	var z7_1_1;
	var x8_1_1;
	var y8_1_1;
	var z8_1_1;
	var x9_1_1;
	var y9_1_1;
	var z9_1_1;
	var x1_2_1;
	var y1_2_1;
	var z1_2_1;
	var x2_2_1;
	var y2_2_1;
	var z2_2_1;
	var x3_2_1;
	var y3_2_1;
	var z3_2_1;
	var x4_2_1;
	var y4_2_1;
	var z4_2_1;
	var x5_2_1;
	var y5_2_1;
	var z5_2_1;
	var x6_2_1;
	var y6_2_1;
	var z6_2_1;
	var x7_2_1;
	var y7_2_1;
	var z7_2_1;
	var x8_2_1;
	var y8_2_1;
	var z8_2_1;
	var x9_2_1;
	var y9_2_1;
	var z9_2_1;
	var x1_3_1;
	var y1_3_1;
	var z1_3_1;
	var x2_3_1;
	var y2_3_1;
	var z2_3_1;
	var x3_3_1;
	var y3_3_1;
	var z3_3_1;
	var x4_3_1;
	var y4_3_1;
	var z4_3_1;
	var x5_3_1;
	var y5_3_1;
	var z5_3_1;
	var x6_3_1;
	var y6_3_1;
	var z6_3_1;
	var x7_3_1;
	var y7_3_1;
	var z7_3_1;
	var x8_3_1;
	var y8_3_1;
	var z8_3_1;
	var x9_3_1;
	var y9_3_1;
	var z9_3_1;
	var x1_4_1;
	var y1_4_1;
	var z1_4_1;
	var x2_4_1;
	var y2_4_1;
	var z2_4_1;
	var x3_4_1;
	var y3_4_1;
	var z3_4_1;
	var x4_4_1;
	var y4_4_1;
	var z4_4_1;
	var x5_4_1;
	var y5_4_1;
	var z5_4_1;
	var x6_4_1;
	var y6_4_1;
	var z6_4_1;
	var x7_4_1;
	var y7_4_1;
	var z7_4_1;
	var x8_4_1;
	var y8_4_1;
	var z8_4_1;
	var x9_4_1;
	var y9_4_1;
	var z9_4_1;
	var x1_5_1;
	var y1_5_1;
	var z1_5_1;
	var x2_5_1;
	var y2_5_1;
	var z2_5_1;
	var x3_5_1;
	var y3_5_1;
	var z3_5_1;
	var x4_5_1;
	var y4_5_1;
	var z4_5_1;
	var x5_5_1;
	var y5_5_1;
	var z5_5_1;
	var x6_5_1;
	var y6_5_1;
	var z6_5_1;
	var x7_5_1;
	var y7_5_1;
	var z7_5_1;
	var x8_5_1;
	var y8_5_1;
	var z8_5_1;
	var x9_5_1;
	var y9_5_1;
	var z9_5_1;
	var x1_6_1;
	var y1_6_1;
	var z1_6_1;
	var x2_6_1;
	var y2_6_1;
	var z2_6_1;
	var x3_6_1;
	var y3_6_1;
	var z3_6_1;
	var x4_6_1;
	var y4_6_1;
	var z4_6_1;
	var x5_6_1;
	var y5_6_1;
	var z5_6_1;
	var x6_6_1;
	var y6_6_1;
	var z6_6_1;
	var x7_6_1;
	var y7_6_1;
	var z7_6_1;
	var x8_6_1;
	var y8_6_1;
	var z8_6_1;
	var x9_6_1;
	var y9_6_1;
	var z9_6_1;
	var x1_7_1;
	var y1_7_1;
	var z1_7_1;
	var x2_7_1;
	var y2_7_1;
	var z2_7_1;
	var x3_7_1;
	var y3_7_1;
	var z3_7_1;
	var x4_7_1;
	var y4_7_1;
	var z4_7_1;
	var x5_7_1;
	var y5_7_1;
	var z5_7_1;
	var x6_7_1;
	var y6_7_1;
	var z6_7_1;
	var x7_7_1;
	var y7_7_1;
	var z7_7_1;
	var x8_7_1;
	var y8_7_1;
	var z8_7_1;
	var x9_7_1;
	var y9_7_1;
	var z9_7_1;
	var x1_8_1;
	var y1_8_1;
	var z1_8_1;
	var x2_8_1;
	var y2_8_1;
	var z2_8_1;
	var x3_8_1;
	var y3_8_1;
	var z3_8_1;
	var x4_8_1;
	var y4_8_1;
	var z4_8_1;
	var x5_8_1;
	var y5_8_1;
	var z5_8_1;
	var x6_8_1;
	var y6_8_1;
	var z6_8_1;
	var x7_8_1;
	var y7_8_1;
	var z7_8_1;
	var x8_8_1;
	var y8_8_1;
	var z8_8_1;
	var x9_8_1;
	var y9_8_1;
	var z9_8_1;
	var x1_9_1;
	var y1_9_1;
	var z1_9_1;
	var x2_9_1;
	var y2_9_1;
	var z2_9_1;
	var x3_9_1;
	var y3_9_1;
	var z3_9_1;
	var x4_9_1;
	var y4_9_1;
	var z4_9_1;
	var x5_9_1;
	var y5_9_1;
	var z5_9_1;
	var x6_9_1;
	var y6_9_1;
	var z6_9_1;
	var x7_9_1;
	var y7_9_1;
	var z7_9_1;
	var x8_9_1;
	var y8_9_1;
	var z8_9_1;
	var x9_9_1;
	var y9_9_1;
	var z9_9_1;
	var x1_1_2;
	var y1_1_2;
	var z1_1_2;
	var x2_1_2;
	var y2_1_2;
	var z2_1_2;
	var x3_1_2;
	var y3_1_2;
	var z3_1_2;
	var x4_1_2;
	var y4_1_2;
	var z4_1_2;
	var x5_1_2;
	var y5_1_2;
	var z5_1_2;
	var x6_1_2;
	var y6_1_2;
	var z6_1_2;
	var x7_1_2;
	var y7_1_2;
	var z7_1_2;
	var x8_1_2;
	var y8_1_2;
	var z8_1_2;
	var x9_1_2;
	var y9_1_2;
	var z9_1_2;
	var x1_2_2;
	var y1_2_2;
	var z1_2_2;
	var x2_2_2;
	var y2_2_2;
	var z2_2_2;
	var x3_2_2;
	var y3_2_2;
	var z3_2_2;
	var x4_2_2;
	var y4_2_2;
	var z4_2_2;
	var x5_2_2;
	var y5_2_2;
	var z5_2_2;
	var x6_2_2;
	var y6_2_2;
	var z6_2_2;
	var x7_2_2;
	var y7_2_2;
	var z7_2_2;
	var x8_2_2;
	var y8_2_2;
	var z8_2_2;
	var x9_2_2;
	var y9_2_2;
	var z9_2_2;
	var x1_3_2;
	var y1_3_2;
	var z1_3_2;
	var x2_3_2;
	var y2_3_2;
	var z2_3_2;
	var x3_3_2;
	var y3_3_2;
	var z3_3_2;
	var x4_3_2;
	var y4_3_2;
	var z4_3_2;
	var x5_3_2;
	var y5_3_2;
	var z5_3_2;
	var x6_3_2;
	var y6_3_2;
	var z6_3_2;
	var x7_3_2;
	var y7_3_2;
	var z7_3_2;
	var x8_3_2;
	var y8_3_2;
	var z8_3_2;
	var x9_3_2;
	var y9_3_2;
	var z9_3_2;
	var x1_4_2;
	var y1_4_2;
	var z1_4_2;
	var x2_4_2;
	var y2_4_2;
	var z2_4_2;
	var x3_4_2;
	var y3_4_2;
	var z3_4_2;
	var x4_4_2;
	var y4_4_2;
	var z4_4_2;
	var x5_4_2;
	var y5_4_2;
	var z5_4_2;
	var x6_4_2;
	var y6_4_2;
	var z6_4_2;
	var x7_4_2;
	var y7_4_2;
	var z7_4_2;
	var x8_4_2;
	var y8_4_2;
	var z8_4_2;
	var x9_4_2;
	var y9_4_2;
	var z9_4_2;
	var x1_5_2;
	var y1_5_2;
	var z1_5_2;
	var x2_5_2;
	var y2_5_2;
	var z2_5_2;
	var x3_5_2;
	var y3_5_2;
	var z3_5_2;
	var x4_5_2;
	var y4_5_2;
	var z4_5_2;
	var x5_5_2;
	var y5_5_2;
	var z5_5_2;
	var x6_5_2;
	var y6_5_2;
	var z6_5_2;
	var x7_5_2;
	var y7_5_2;
	var z7_5_2;
	var x8_5_2;
	var y8_5_2;
	var z8_5_2;
	var x9_5_2;
	var y9_5_2;
	var z9_5_2;
	var x1_6_2;
	var y1_6_2;
	var z1_6_2;
	var x2_6_2;
	var y2_6_2;
	var z2_6_2;
	var x3_6_2;
	var y3_6_2;
	var z3_6_2;
	var x4_6_2;
	var y4_6_2;
	var z4_6_2;
	var x5_6_2;
	var y5_6_2;
	var z5_6_2;
	var x6_6_2;
	var y6_6_2;
	var z6_6_2;
	var x7_6_2;
	var y7_6_2;
	var z7_6_2;
	var x8_6_2;
	var y8_6_2;
	var z8_6_2;
	var x9_6_2;
	var y9_6_2;
	var z9_6_2;
	var x1_7_2;
	var y1_7_2;
	var z1_7_2;
	var x2_7_2;
	var y2_7_2;
	var z2_7_2;
	var x3_7_2;
	var y3_7_2;
	var z3_7_2;
	var x4_7_2;
	var y4_7_2;
	var z4_7_2;
	var x5_7_2;
	var y5_7_2;
	var z5_7_2;
	var x6_7_2;
	var y6_7_2;
	var z6_7_2;
	var x7_7_2;
	var y7_7_2;
	var z7_7_2;
	var x8_7_2;
	var y8_7_2;
	var z8_7_2;
	var x9_7_2;
	var y9_7_2;
	var z9_7_2;
	var x1_8_2;
	var y1_8_2;
	var z1_8_2;
	var x2_8_2;
	var y2_8_2;
	var z2_8_2;
	var x3_8_2;
	var y3_8_2;
	var z3_8_2;
	var x4_8_2;
	var y4_8_2;
	var z4_8_2;
	var x5_8_2;
	var y5_8_2;
	var z5_8_2;
	var x6_8_2;
	var y6_8_2;
	var z6_8_2;
	var x7_8_2;
	var y7_8_2;
	var z7_8_2;
	var x8_8_2;
	var y8_8_2;
	var z8_8_2;
	var x9_8_2;
	var y9_8_2;
	var z9_8_2;
	var x1_9_2;
	var y1_9_2;
	var z1_9_2;
	var x2_9_2;
	var y2_9_2;
	var z2_9_2;
	var x3_9_2;
	var y3_9_2;
	var z3_9_2;
	var x4_9_2;
	var y4_9_2;
	var z4_9_2;
	var x5_9_2;
	var y5_9_2;
	var z5_9_2;
	var x6_9_2;
	var y6_9_2;
	var z6_9_2;
	var x7_9_2;
	var y7_9_2;
	var z7_9_2;
	var x8_9_2;
	var y8_9_2;
	var z8_9_2;
	var x9_9_2;
	var y9_9_2;
	var z9_9_2;
	var x1_1_3;
	var y1_1_3;
	var z1_1_3;
	var x2_1_3;
	var y2_1_3;
	var z2_1_3;
	var x3_1_3;
	var y3_1_3;
	var z3_1_3;
	var x4_1_3;
	var y4_1_3;
	var z4_1_3;
	var x5_1_3;
	var y5_1_3;
	var z5_1_3;
	var x6_1_3;
	var y6_1_3;
	var z6_1_3;
	var x7_1_3;
	var y7_1_3;
	var z7_1_3;
	var x8_1_3;
	var y8_1_3;
	var z8_1_3;
	var x9_1_3;
	var y9_1_3;
	var z9_1_3;
	var x1_2_3;
	var y1_2_3;
	var z1_2_3;
	var x2_2_3;
	var y2_2_3;
	var z2_2_3;
	var x3_2_3;
	var y3_2_3;
	var z3_2_3;
	var x4_2_3;
	var y4_2_3;
	var z4_2_3;
	var x5_2_3;
	var y5_2_3;
	var z5_2_3;
	var x6_2_3;
	var y6_2_3;
	var z6_2_3;
	var x7_2_3;
	var y7_2_3;
	var z7_2_3;
	var x8_2_3;
	var y8_2_3;
	var z8_2_3;
	var x9_2_3;
	var y9_2_3;
	var z9_2_3;
	var x1_3_3;
	var y1_3_3;
	var z1_3_3;
	var x2_3_3;
	var y2_3_3;
	var z2_3_3;
	var x3_3_3;
	var y3_3_3;
	var z3_3_3;
	var x4_3_3;
	var y4_3_3;
	var z4_3_3;
	var x5_3_3;
	var y5_3_3;
	var z5_3_3;
	var x6_3_3;
	var y6_3_3;
	var z6_3_3;
	var x7_3_3;
	var y7_3_3;
	var z7_3_3;
	var x8_3_3;
	var y8_3_3;
	var z8_3_3;
	var x9_3_3;
	var y9_3_3;
	var z9_3_3;
	var x1_4_3;
	var y1_4_3;
	var z1_4_3;
	var x2_4_3;
	var y2_4_3;
	var z2_4_3;
	var x3_4_3;
	var y3_4_3;
	var z3_4_3;
	var x4_4_3;
	var y4_4_3;
	var z4_4_3;
	var x5_4_3;
	var y5_4_3;
	var z5_4_3;
	var x6_4_3;
	var y6_4_3;
	var z6_4_3;
	var x7_4_3;
	var y7_4_3;
	var z7_4_3;
	var x8_4_3;
	var y8_4_3;
	var z8_4_3;
	var x9_4_3;
	var y9_4_3;
	var z9_4_3;
	var x1_5_3;
	var y1_5_3;
	var z1_5_3;
	var x2_5_3;
	var y2_5_3;
	var z2_5_3;
	var x3_5_3;
	var y3_5_3;
	var z3_5_3;
	var x4_5_3;
	var y4_5_3;
	var z4_5_3;
	var x5_5_3;
	var y5_5_3;
	var z5_5_3;
	var x6_5_3;
	var y6_5_3;
	var z6_5_3;
	var x7_5_3;
	var y7_5_3;
	var z7_5_3;
	var x8_5_3;
	var y8_5_3;
	var z8_5_3;
	var x9_5_3;
	var y9_5_3;
	var z9_5_3;
	var x1_6_3;
	var y1_6_3;
	var z1_6_3;
	var x2_6_3;
	var y2_6_3;
	var z2_6_3;
	var x3_6_3;
	var y3_6_3;
	var z3_6_3;
	var x4_6_3;
	var y4_6_3;
	var z4_6_3;
	var x5_6_3;
	var y5_6_3;
	var z5_6_3;
	var x6_6_3;
	var y6_6_3;
	var z6_6_3;
	var x7_6_3;
	var y7_6_3;
	var z7_6_3;
	var x8_6_3;
	var y8_6_3;
	var z8_6_3;
	var x9_6_3;
	var y9_6_3;
	var z9_6_3;
	var x1_7_3;
	var y1_7_3;
	var z1_7_3;
	var x2_7_3;
	var y2_7_3;
	var z2_7_3;
	var x3_7_3;
	var y3_7_3;
	var z3_7_3;
	var x4_7_3;
	var y4_7_3;
	var z4_7_3;
	var x5_7_3;
	var y5_7_3;
	var z5_7_3;
	var x6_7_3;
	var y6_7_3;
	var z6_7_3;
	var x7_7_3;
	var y7_7_3;
	var z7_7_3;
	var x8_7_3;
	var y8_7_3;
	var z8_7_3;
	var x9_7_3;
	var y9_7_3;
	var z9_7_3;
	var x1_8_3;
	var y1_8_3;
	var z1_8_3;
	var x2_8_3;
	var y2_8_3;
	var z2_8_3;
	var x3_8_3;
	var y3_8_3;
	var z3_8_3;
	var x4_8_3;
	var y4_8_3;
	var z4_8_3;
	var x5_8_3;
	var y5_8_3;
	var z5_8_3;
	var x6_8_3;
	var y6_8_3;
	var z6_8_3;
	var x7_8_3;
	var y7_8_3;
	var z7_8_3;
	var x8_8_3;
	var y8_8_3;
	var z8_8_3;
	var x9_8_3;
	var y9_8_3;
	var z9_8_3;
	var x1_9_3;
	var y1_9_3;
	var z1_9_3;
	var x2_9_3;
	var y2_9_3;
	var z2_9_3;
	var x3_9_3;
	var y3_9_3;
	var z3_9_3;
	var x4_9_3;
	var y4_9_3;
	var z4_9_3;
	var x5_9_3;
	var y5_9_3;
	var z5_9_3;
	var x6_9_3;
	var y6_9_3;
	var z6_9_3;
	var x7_9_3;
	var y7_9_3;
	var z7_9_3;
	var x8_9_3;
	var y8_9_3;
	var z8_9_3;
	var x9_9_3;
	var y9_9_3;
	var z9_9_3;
	var x1_1_4;
	var y1_1_4;
	var z1_1_4;
	var x2_1_4;
	var y2_1_4;
	var z2_1_4;
	var x3_1_4;
	var y3_1_4;
	var z3_1_4;
	var x4_1_4;
	var y4_1_4;
	var z4_1_4;
	var x5_1_4;
	var y5_1_4;
	var z5_1_4;
	var x6_1_4;
	var y6_1_4;
	var z6_1_4;
	var x7_1_4;
	var y7_1_4;
	var z7_1_4;
	var x8_1_4;
	var y8_1_4;
	var z8_1_4;
	var x9_1_4;
	var y9_1_4;
	var z9_1_4;
	var x1_2_4;
	var y1_2_4;
	var z1_2_4;
	var x2_2_4;
	var y2_2_4;
	var z2_2_4;
	var x3_2_4;
	var y3_2_4;
	var z3_2_4;
	var x4_2_4;
	var y4_2_4;
	var z4_2_4;
	var x5_2_4;
	var y5_2_4;
	var z5_2_4;
	var x6_2_4;
	var y6_2_4;
	var z6_2_4;
	var x7_2_4;
	var y7_2_4;
	var z7_2_4;
	var x8_2_4;
	var y8_2_4;
	var z8_2_4;
	var x9_2_4;
	var y9_2_4;
	var z9_2_4;
	var x1_3_4;
	var y1_3_4;
	var z1_3_4;
	var x2_3_4;
	var y2_3_4;
	var z2_3_4;
	var x3_3_4;
	var y3_3_4;
	var z3_3_4;
	var x4_3_4;
	var y4_3_4;
	var z4_3_4;
	var x5_3_4;
	var y5_3_4;
	var z5_3_4;
	var x6_3_4;
	var y6_3_4;
	var z6_3_4;
	var x7_3_4;
	var y7_3_4;
	var z7_3_4;
	var x8_3_4;
	var y8_3_4;
	var z8_3_4;
	var x9_3_4;
	var y9_3_4;
	var z9_3_4;
	var x1_4_4;
	var y1_4_4;
	var z1_4_4;
	var x2_4_4;
	var y2_4_4;
	var z2_4_4;
	var x3_4_4;
	var y3_4_4;
	var z3_4_4;
	var x4_4_4;
	var y4_4_4;
	var z4_4_4;
	var x5_4_4;
	var y5_4_4;
	var z5_4_4;
	var x6_4_4;
	var y6_4_4;
	var z6_4_4;
	var x7_4_4;
	var y7_4_4;
	var z7_4_4;
	var x8_4_4;
	var y8_4_4;
	var z8_4_4;
	var x9_4_4;
	var y9_4_4;
	var z9_4_4;
	var x1_5_4;
	var y1_5_4;
	var z1_5_4;
	var x2_5_4;
	var y2_5_4;
	var z2_5_4;
	var x3_5_4;
	var y3_5_4;
	var z3_5_4;
	var x4_5_4;
	var y4_5_4;
	var z4_5_4;
	var x5_5_4;
	var y5_5_4;
	var z5_5_4;
	var x6_5_4;
	var y6_5_4;
	var z6_5_4;
	var x7_5_4;
	var y7_5_4;
	var z7_5_4;
	var x8_5_4;
	var y8_5_4;
	var z8_5_4;
	var x9_5_4;
	var y9_5_4;
	var z9_5_4;
	var x1_6_4;
	var y1_6_4;
	var z1_6_4;
	var x2_6_4;
	var y2_6_4;
	var z2_6_4;
	var x3_6_4;
	var y3_6_4;
	var z3_6_4;
	var x4_6_4;
	var y4_6_4;
	var z4_6_4;
	var x5_6_4;
	var y5_6_4;
	var z5_6_4;
	var x6_6_4;
	var y6_6_4;
	var z6_6_4;
	var x7_6_4;
	var y7_6_4;
	var z7_6_4;
	var x8_6_4;
	var y8_6_4;
	var z8_6_4;
	var x9_6_4;
	var y9_6_4;
	var z9_6_4;
	var x1_7_4;
	var y1_7_4;
	var z1_7_4;
	var x2_7_4;
	var y2_7_4;
	var z2_7_4;
	var x3_7_4;
	var y3_7_4;
	var z3_7_4;
	var x4_7_4;
	var y4_7_4;
	var z4_7_4;
	var x5_7_4;
	var y5_7_4;
	var z5_7_4;
	var x6_7_4;
	var y6_7_4;
	var z6_7_4;
	var x7_7_4;
	var y7_7_4;
	var z7_7_4;
	var x8_7_4;
	var y8_7_4;
	var z8_7_4;
	var x9_7_4;
	var y9_7_4;
	var z9_7_4;
	var x1_8_4;
	var y1_8_4;
	var z1_8_4;
	var x2_8_4;
	var y2_8_4;
	var z2_8_4;
	var x3_8_4;
	var y3_8_4;
	var z3_8_4;
	var x4_8_4;
	var y4_8_4;
	var z4_8_4;
	var x5_8_4;
	var y5_8_4;
	var z5_8_4;
	var x6_8_4;
	var y6_8_4;
	var z6_8_4;
	var x7_8_4;
	var y7_8_4;
	var z7_8_4;
	var x8_8_4;
	var y8_8_4;
	var z8_8_4;
	var x9_8_4;
	var y9_8_4;
	var z9_8_4;
	var x1_9_4;
	var y1_9_4;
	var z1_9_4;
	var x2_9_4;
	var y2_9_4;
	var z2_9_4;
	var x3_9_4;
	var y3_9_4;
	var z3_9_4;
	var x4_9_4;
	var y4_9_4;
	var z4_9_4;
	var x5_9_4;
	var y5_9_4;
	var z5_9_4;
	var x6_9_4;
	var y6_9_4;
	var z6_9_4;
	var x7_9_4;
	var y7_9_4;
	var z7_9_4;
	var x8_9_4;
	var y8_9_4;
	var z8_9_4;
	var x9_9_4;
	var y9_9_4;
	var z9_9_4;
	var x1_1_5;
	var y1_1_5;
	var z1_1_5;
	var x2_1_5;
	var y2_1_5;
	var z2_1_5;
	var x3_1_5;
	var y3_1_5;
	var z3_1_5;
	var x4_1_5;
	var y4_1_5;
	var z4_1_5;
	var x5_1_5;
	var y5_1_5;
	var z5_1_5;
	var x6_1_5;
	var y6_1_5;
	var z6_1_5;
	var x7_1_5;
	var y7_1_5;
	var z7_1_5;
	var x8_1_5;
	var y8_1_5;
	var z8_1_5;
	var x9_1_5;
	var y9_1_5;
	var z9_1_5;
	var x1_2_5;
	var y1_2_5;
	var z1_2_5;
	var x2_2_5;
	var y2_2_5;
	var z2_2_5;
	var x3_2_5;
	var y3_2_5;
	var z3_2_5;
	var x4_2_5;
	var y4_2_5;
	var z4_2_5;
	var x5_2_5;
	var y5_2_5;
	var z5_2_5;
	var x6_2_5;
	var y6_2_5;
	var z6_2_5;
	var x7_2_5;
	var y7_2_5;
	var z7_2_5;
	var x8_2_5;
	var y8_2_5;
	var z8_2_5;
	var x9_2_5;
	var y9_2_5;
	var z9_2_5;
	var x1_3_5;
	var y1_3_5;
	var z1_3_5;
	var x2_3_5;
	var y2_3_5;
	var z2_3_5;
	var x3_3_5;
	var y3_3_5;
	var z3_3_5;
	var x4_3_5;
	var y4_3_5;
	var z4_3_5;
	var x5_3_5;
	var y5_3_5;
	var z5_3_5;
	var x6_3_5;
	var y6_3_5;
	var z6_3_5;
	var x7_3_5;
	var y7_3_5;
	var z7_3_5;
	var x8_3_5;
	var y8_3_5;
	var z8_3_5;
	var x9_3_5;
	var y9_3_5;
	var z9_3_5;
	var x1_4_5;
	var y1_4_5;
	var z1_4_5;
	var x2_4_5;
	var y2_4_5;
	var z2_4_5;
	var x3_4_5;
	var y3_4_5;
	var z3_4_5;
	var x4_4_5;
	var y4_4_5;
	var z4_4_5;
	var x5_4_5;
	var y5_4_5;
	var z5_4_5;
	var x6_4_5;
	var y6_4_5;
	var z6_4_5;
	var x7_4_5;
	var y7_4_5;
	var z7_4_5;
	var x8_4_5;
	var y8_4_5;
	var z8_4_5;
	var x9_4_5;
	var y9_4_5;
	var z9_4_5;
	var x1_5_5;
	var y1_5_5;
	var z1_5_5;
	var x2_5_5;
	var y2_5_5;
	var z2_5_5;
	var x3_5_5;
	var y3_5_5;
	var z3_5_5;
	var x4_5_5;
	var y4_5_5;
	var z4_5_5;
	var x5_5_5;
	var y5_5_5;
	var z5_5_5;
	var x6_5_5;
	var y6_5_5;
	var z6_5_5;
	var x7_5_5;
	var y7_5_5;
	var z7_5_5;
	var x8_5_5;
	var y8_5_5;
	var z8_5_5;
	var x9_5_5;
	var y9_5_5;
	var z9_5_5;
	var x1_6_5;
	var y1_6_5;
	var z1_6_5;
	var x2_6_5;
	var y2_6_5;
	var z2_6_5;
	var x3_6_5;
	var y3_6_5;
	var z3_6_5;
	var x4_6_5;
	var y4_6_5;
	var z4_6_5;
	var x5_6_5;
	var y5_6_5;
	var z5_6_5;
	var x6_6_5;
	var y6_6_5;
	var z6_6_5;
	var x7_6_5;
	var y7_6_5;
	var z7_6_5;
	var x8_6_5;
	var y8_6_5;
	var z8_6_5;
	var x9_6_5;
	var y9_6_5;
	var z9_6_5;
	var x1_7_5;
	var y1_7_5;
	var z1_7_5;
	var x2_7_5;
	var y2_7_5;
	var z2_7_5;
	var x3_7_5;
	var y3_7_5;
	var z3_7_5;
	var x4_7_5;
	var y4_7_5;
	var z4_7_5;
	var x5_7_5;
	var y5_7_5;
	var z5_7_5;
	var x6_7_5;
	var y6_7_5;
	var z6_7_5;
	var x7_7_5;
	var y7_7_5;
	var z7_7_5;
	var x8_7_5;
	var y8_7_5;
	var z8_7_5;
	var x9_7_5;
	var y9_7_5;
	var z9_7_5;
	var x1_8_5;
	var y1_8_5;
	var z1_8_5;
	var x2_8_5;
	var y2_8_5;
	var z2_8_5;
	var x3_8_5;
	var y3_8_5;
	var z3_8_5;
	var x4_8_5;
	var y4_8_5;
	var z4_8_5;
	var x5_8_5;
	var y5_8_5;
	var z5_8_5;
	var x6_8_5;
	var y6_8_5;
	var z6_8_5;
	var x7_8_5;
	var y7_8_5;
	var z7_8_5;
	var x8_8_5;
	var y8_8_5;
	var z8_8_5;
	var x9_8_5;
	var y9_8_5;
	var z9_8_5;
	var x1_9_5;
	var y1_9_5;
	var z1_9_5;
	var x2_9_5;
	var y2_9_5;
	var z2_9_5;
	var x3_9_5;
	var y3_9_5;
	var z3_9_5;
	var x4_9_5;
	var y4_9_5;
	var z4_9_5;
	var x5_9_5;
	var y5_9_5;
	var z5_9_5;
	var x6_9_5;
	var y6_9_5;
	var z6_9_5;
	var x7_9_5;
	var y7_9_5;
	var z7_9_5;
	var x8_9_5;
	var y8_9_5;
	var z8_9_5;
	var x9_9_5;
	var y9_9_5;
	var z9_9_5;
	var x1_1_6;
	var y1_1_6;
	var z1_1_6;
	var x2_1_6;
	var y2_1_6;
	var z2_1_6;
	var x3_1_6;
	var y3_1_6;
	var z3_1_6;
	var x4_1_6;
	var y4_1_6;
	var z4_1_6;
	var x5_1_6;
	var y5_1_6;
	var z5_1_6;
	var x6_1_6;
	var y6_1_6;
	var z6_1_6;
	var x7_1_6;
	var y7_1_6;
	var z7_1_6;
	var x8_1_6;
	var y8_1_6;
	var z8_1_6;
	var x9_1_6;
	var y9_1_6;
	var z9_1_6;
	var x1_2_6;
	var y1_2_6;
	var z1_2_6;
	var x2_2_6;
	var y2_2_6;
	var z2_2_6;
	var x3_2_6;
	var y3_2_6;
	var z3_2_6;
	var x4_2_6;
	var y4_2_6;
	var z4_2_6;
	var x5_2_6;
	var y5_2_6;
	var z5_2_6;
	var x6_2_6;
	var y6_2_6;
	var z6_2_6;
	var x7_2_6;
	var y7_2_6;
	var z7_2_6;
	var x8_2_6;
	var y8_2_6;
	var z8_2_6;
	var x9_2_6;
	var y9_2_6;
	var z9_2_6;
	var x1_3_6;
	var y1_3_6;
	var z1_3_6;
	var x2_3_6;
	var y2_3_6;
	var z2_3_6;
	var x3_3_6;
	var y3_3_6;
	var z3_3_6;
	var x4_3_6;
	var y4_3_6;
	var z4_3_6;
	var x5_3_6;
	var y5_3_6;
	var z5_3_6;
	var x6_3_6;
	var y6_3_6;
	var z6_3_6;
	var x7_3_6;
	var y7_3_6;
	var z7_3_6;
	var x8_3_6;
	var y8_3_6;
	var z8_3_6;
	var x9_3_6;
	var y9_3_6;
	var z9_3_6;
	var x1_4_6;
	var y1_4_6;
	var z1_4_6;
	var x2_4_6;
	var y2_4_6;
	var z2_4_6;
	var x3_4_6;
	var y3_4_6;
	var z3_4_6;
	var x4_4_6;
	var y4_4_6;
	var z4_4_6;
	var x5_4_6;
	var y5_4_6;
	var z5_4_6;
	var x6_4_6;
	var y6_4_6;
	var z6_4_6;
	var x7_4_6;
	var y7_4_6;
	var z7_4_6;
	var x8_4_6;
	var y8_4_6;
	var z8_4_6;
	var x9_4_6;
	var y9_4_6;
	var z9_4_6;
	var x1_5_6;
	var y1_5_6;
	var z1_5_6;
	var x2_5_6;
	var y2_5_6;
	var z2_5_6;
	var x3_5_6;
	var y3_5_6;
	var z3_5_6;
	var x4_5_6;
	var y4_5_6;
	var z4_5_6;
	var x5_5_6;
	var y5_5_6;
	var z5_5_6;
	var x6_5_6;
	var y6_5_6;
	var z6_5_6;
	var x7_5_6;
	var y7_5_6;
	var z7_5_6;
	var x8_5_6;
	var y8_5_6;
	var z8_5_6;
	var x9_5_6;
	var y9_5_6;
	var z9_5_6;
	var x1_6_6;
	var y1_6_6;
	var z1_6_6;
	var x2_6_6;
	var y2_6_6;
	var z2_6_6;
	var x3_6_6;
	var y3_6_6;
	var z3_6_6;
	var x4_6_6;
	var y4_6_6;
	var z4_6_6;
	var x5_6_6;
	var y5_6_6;
	var z5_6_6;
	var x6_6_6;
	var y6_6_6;
	var z6_6_6;
	var x7_6_6;
	var y7_6_6;
	var z7_6_6;
	var x8_6_6;
	var y8_6_6;
	var z8_6_6;
	var x9_6_6;
	var y9_6_6;
	var z9_6_6;
	var x1_7_6;
	var y1_7_6;
	var z1_7_6;
	var x2_7_6;
	var y2_7_6;
	var z2_7_6;
	var x3_7_6;
	var y3_7_6;
	var z3_7_6;
	var x4_7_6;
	var y4_7_6;
	var z4_7_6;
	var x5_7_6;
	var y5_7_6;
	var z5_7_6;
	var x6_7_6;
	var y6_7_6;
	var z6_7_6;
	var x7_7_6;
	var y7_7_6;
	var z7_7_6;
	var x8_7_6;
	var y8_7_6;
	var z8_7_6;
	var x9_7_6;
	var y9_7_6;
	var z9_7_6;
	var x1_8_6;
	var y1_8_6;
	var z1_8_6;
	var x2_8_6;
	var y2_8_6;
	var z2_8_6;
	var x3_8_6;
	var y3_8_6;
	var z3_8_6;
	var x4_8_6;
	var y4_8_6;
	var z4_8_6;
	var x5_8_6;
	var y5_8_6;
	var z5_8_6;
	var x6_8_6;
	var y6_8_6;
	var z6_8_6;
	var x7_8_6;
	var y7_8_6;
	var z7_8_6;
	var x8_8_6;
	var y8_8_6;
	var z8_8_6;
	var x9_8_6;
	var y9_8_6;
	var z9_8_6;
	var x1_9_6;
	var y1_9_6;
	var z1_9_6;
	var x2_9_6;
	var y2_9_6;
	var z2_9_6;
	var x3_9_6;
	var y3_9_6;
	var z3_9_6;
	var x4_9_6;
	var y4_9_6;
	var z4_9_6;
	var x5_9_6;
	var y5_9_6;
	var z5_9_6;
	var x6_9_6;
	var y6_9_6;
	var z6_9_6;
	var x7_9_6;
	var y7_9_6;
	var z7_9_6;
	var x8_9_6;
	var y8_9_6;
	var z8_9_6;
	var x9_9_6;
	var y9_9_6;
	var z9_9_6;
	var x1_1_7;
	var y1_1_7;
	var z1_1_7;
	var x2_1_7;
	var y2_1_7;
	var z2_1_7;
	var x3_1_7;
	var y3_1_7;
	var z3_1_7;
	var x4_1_7;
	var y4_1_7;
	var z4_1_7;
	var x5_1_7;
	var y5_1_7;
	var z5_1_7;
	var x6_1_7;
	var y6_1_7;
	var z6_1_7;
	var x7_1_7;
	var y7_1_7;
	var z7_1_7;
	var x8_1_7;
	var y8_1_7;
	var z8_1_7;
	var x9_1_7;
	var y9_1_7;
	var z9_1_7;
	var x1_2_7;
	var y1_2_7;
	var z1_2_7;
	var x2_2_7;
	var y2_2_7;
	var z2_2_7;
	var x3_2_7;
	var y3_2_7;
	var z3_2_7;
	var x4_2_7;
	var y4_2_7;
	var z4_2_7;
	var x5_2_7;
	var y5_2_7;
	var z5_2_7;
	var x6_2_7;
	var y6_2_7;
	var z6_2_7;
	var x7_2_7;
	var y7_2_7;
	var z7_2_7;
	var x8_2_7;
	var y8_2_7;
	var z8_2_7;
	var x9_2_7;
	var y9_2_7;
	var z9_2_7;
	var x1_3_7;
	var y1_3_7;
	var z1_3_7;
	var x2_3_7;
	var y2_3_7;
	var z2_3_7;
	var x3_3_7;
	var y3_3_7;
	var z3_3_7;
	var x4_3_7;
	var y4_3_7;
	var z4_3_7;
	var x5_3_7;
	var y5_3_7;
	var z5_3_7;
	var x6_3_7;
	var y6_3_7;
	var z6_3_7;
	var x7_3_7;
	var y7_3_7;
	var z7_3_7;
	var x8_3_7;
	var y8_3_7;
	var z8_3_7;
	var x9_3_7;
	var y9_3_7;
	var z9_3_7;
	var x1_4_7;
	var y1_4_7;
	var z1_4_7;
	var x2_4_7;
	var y2_4_7;
	var z2_4_7;
	var x3_4_7;
	var y3_4_7;
	var z3_4_7;
	var x4_4_7;
	var y4_4_7;
	var z4_4_7;
	var x5_4_7;
	var y5_4_7;
	var z5_4_7;
	var x6_4_7;
	var y6_4_7;
	var z6_4_7;
	var x7_4_7;
	var y7_4_7;
	var z7_4_7;
	var x8_4_7;
	var y8_4_7;
	var z8_4_7;
	var x9_4_7;
	var y9_4_7;
	var z9_4_7;
	var x1_5_7;
	var y1_5_7;
	var z1_5_7;
	var x2_5_7;
	var y2_5_7;
	var z2_5_7;
	var x3_5_7;
	var y3_5_7;
	var z3_5_7;
	var x4_5_7;
	var y4_5_7;
	var z4_5_7;
	var x5_5_7;
	var y5_5_7;
	var z5_5_7;
	var x6_5_7;
	var y6_5_7;
	var z6_5_7;
	var x7_5_7;
	var y7_5_7;
	var z7_5_7;
	var x8_5_7;
	var y8_5_7;
	var z8_5_7;
	var x9_5_7;
	var y9_5_7;
	var z9_5_7;
	var x1_6_7;
	var y1_6_7;
	var z1_6_7;
	var x2_6_7;
	var y2_6_7;
	var z2_6_7;
	var x3_6_7;
	var y3_6_7;
	var z3_6_7;
	var x4_6_7;
	var y4_6_7;
	var z4_6_7;
	var x5_6_7;
	var y5_6_7;
	var z5_6_7;
	var x6_6_7;
	var y6_6_7;
	var z6_6_7;
	var x7_6_7;
	var y7_6_7;
	var z7_6_7;
	var x8_6_7;
	var y8_6_7;
	var z8_6_7;
	var x9_6_7;
	var y9_6_7;
	var z9_6_7;
	var x1_7_7;
	var y1_7_7;
	var z1_7_7;
	var x2_7_7;
	var y2_7_7;
	var z2_7_7;
	var x3_7_7;
	var y3_7_7;
	var z3_7_7;
	var x4_7_7;
	var y4_7_7;
	var z4_7_7;
	var x5_7_7;
	var y5_7_7;
	var z5_7_7;
	var x6_7_7;
	var y6_7_7;
	var z6_7_7;
	var x7_7_7;
	var y7_7_7;
	var z7_7_7;
	var x8_7_7;
	var y8_7_7;
	var z8_7_7;
	var x9_7_7;
	var y9_7_7;
	var z9_7_7;
	var x1_8_7;
	var y1_8_7;
	var z1_8_7;
	var x2_8_7;
	var y2_8_7;
	var z2_8_7;
	var x3_8_7;
	var y3_8_7;
	var z3_8_7;
	var x4_8_7;
	var y4_8_7;
	var z4_8_7;
	var x5_8_7;
	var y5_8_7;
	var z5_8_7;
	var x6_8_7;
	var y6_8_7;
	var z6_8_7;
	var x7_8_7;
	var y7_8_7;
	var z7_8_7;
	var x8_8_7;
	var y8_8_7;
	var z8_8_7;
	var x9_8_7;
	var y9_8_7;
	var z9_8_7;
	var x1_9_7;
	var y1_9_7;
	var z1_9_7;
	var x2_9_7;
	var y2_9_7;
	var z2_9_7;
	var x3_9_7;
	var y3_9_7;
	var z3_9_7;
	var x4_9_7;
	var y4_9_7;
	var z4_9_7;
	var x5_9_7;
	var y5_9_7;
	var z5_9_7;
	var x6_9_7;
	var y6_9_7;
	var z6_9_7;
	var x7_9_7;
	var y7_9_7;
	var z7_9_7;
	var x8_9_7;
	var y8_9_7;
	var z8_9_7;
	var x9_9_7;
	var y9_9_7;
	var z9_9_7;
	var x1_1_8;
	var y1_1_8;
	var z1_1_8;
	var x2_1_8;
	var y2_1_8;
	var z2_1_8;
	var x3_1_8;
	var y3_1_8;
	var z3_1_8;
	var x4_1_8;
	var y4_1_8;
	var z4_1_8;
	var x5_1_8;
	var y5_1_8;
	var z5_1_8;
	var x6_1_8;
	var y6_1_8;
	var z6_1_8;
	var x7_1_8;
	var y7_1_8;
	var z7_1_8;
	var x8_1_8;
	var y8_1_8;
	var z8_1_8;
	var x9_1_8;
	var y9_1_8;
	var z9_1_8;
	var x1_2_8;
	var y1_2_8;
	var z1_2_8;
	var x2_2_8;
	var y2_2_8;
	var z2_2_8;
	var x3_2_8;
	var y3_2_8;
	var z3_2_8;
	var x4_2_8;
	var y4_2_8;
	var z4_2_8;
	var x5_2_8;
	var y5_2_8;
	var z5_2_8;
	var x6_2_8;
	var y6_2_8;
	var z6_2_8;
	var x7_2_8;
	var y7_2_8;
	var z7_2_8;
	var x8_2_8;
	var y8_2_8;
	var z8_2_8;
	var x9_2_8;
	var y9_2_8;
	var z9_2_8;
	var x1_3_8;
	var y1_3_8;
	var z1_3_8;
	var x2_3_8;
	var y2_3_8;
	var z2_3_8;
	var x3_3_8;
	var y3_3_8;
	var z3_3_8;
	var x4_3_8;
	var y4_3_8;
	var z4_3_8;
	var x5_3_8;
	var y5_3_8;
	var z5_3_8;
	var x6_3_8;
	var y6_3_8;
	var z6_3_8;
	var x7_3_8;
	var y7_3_8;
	var z7_3_8;
	var x8_3_8;
	var y8_3_8;
	var z8_3_8;
	var x9_3_8;
	var y9_3_8;
	var z9_3_8;
	var x1_4_8;
	var y1_4_8;
	var z1_4_8;
	var x2_4_8;
	var y2_4_8;
	var z2_4_8;
	var x3_4_8;
	var y3_4_8;
	var z3_4_8;
	var x4_4_8;
	var y4_4_8;
	var z4_4_8;
	var x5_4_8;
	var y5_4_8;
	var z5_4_8;
	var x6_4_8;
	var y6_4_8;
	var z6_4_8;
	var x7_4_8;
	var y7_4_8;
	var z7_4_8;
	var x8_4_8;
	var y8_4_8;
	var z8_4_8;
	var x9_4_8;
	var y9_4_8;
	var z9_4_8;
	var x1_5_8;
	var y1_5_8;
	var z1_5_8;
	var x2_5_8;
	var y2_5_8;
	var z2_5_8;
	var x3_5_8;
	var y3_5_8;
	var z3_5_8;
	var x4_5_8;
	var y4_5_8;
	var z4_5_8;
	var x5_5_8;
	var y5_5_8;
	var z5_5_8;
	var x6_5_8;
	var y6_5_8;
	var z6_5_8;
	var x7_5_8;
	var y7_5_8;
	var z7_5_8;
	var x8_5_8;
	var y8_5_8;
	var z8_5_8;
	var x9_5_8;
	var y9_5_8;
	var z9_5_8;
	var x1_6_8;
	var y1_6_8;
	var z1_6_8;
	var x2_6_8;
	var y2_6_8;
	var z2_6_8;
	var x3_6_8;
	var y3_6_8;
	var z3_6_8;
	var x4_6_8;
	var y4_6_8;
	var z4_6_8;
	var x5_6_8;
	var y5_6_8;
	var z5_6_8;
	var x6_6_8;
	var y6_6_8;
	var z6_6_8;
	var x7_6_8;
	var y7_6_8;
	var z7_6_8;
	var x8_6_8;
	var y8_6_8;
	var z8_6_8;
	var x9_6_8;
	var y9_6_8;
	var z9_6_8;
	var x1_7_8;
	var y1_7_8;
	var z1_7_8;
	var x2_7_8;
	var y2_7_8;
	var z2_7_8;
	var x3_7_8;
	var y3_7_8;
	var z3_7_8;
	var x4_7_8;
	var y4_7_8;
	var z4_7_8;
	var x5_7_8;
	var y5_7_8;
	var z5_7_8;
	var x6_7_8;
	var y6_7_8;
	var z6_7_8;
	var x7_7_8;
	var y7_7_8;
	var z7_7_8;
	var x8_7_8;
	var y8_7_8;
	var z8_7_8;
	var x9_7_8;
	var y9_7_8;
	var z9_7_8;
	var x1_8_8;
	var y1_8_8;
	var z1_8_8;
	var x2_8_8;
	var y2_8_8;
	var z2_8_8;
	var x3_8_8;
	var y3_8_8;
	var z3_8_8;
	var x4_8_8;
	var y4_8_8;
	var z4_8_8;
	var x5_8_8;
	var y5_8_8;
	var z5_8_8;
	var x6_8_8;
	var y6_8_8;
	var z6_8_8;
	var x7_8_8;
	var y7_8_8;
	var z7_8_8;
	var x8_8_8;
	var y8_8_8;
	var z8_8_8;
	var x9_8_8;
	var y9_8_8;
	var z9_8_8;
	var x1_9_8;
	var y1_9_8;
	var z1_9_8;
	var x2_9_8;
	var y2_9_8;
	var z2_9_8;
	var x3_9_8;
	var y3_9_8;
	var z3_9_8;
	var x4_9_8;
	var y4_9_8;
	var z4_9_8;
	var x5_9_8;
	var y5_9_8;
	var z5_9_8;
	var x6_9_8;
	var y6_9_8;
	var z6_9_8;
	var x7_9_8;
	var y7_9_8;
	var z7_9_8;
	var x8_9_8;
	var y8_9_8;
	var z8_9_8;
	var x9_9_8;
	var y9_9_8;
	var z9_9_8;
	var x1_1_9;
	var y1_1_9;
	var z1_1_9;
	var x2_1_9;
	var y2_1_9;
	var z2_1_9;
	var x3_1_9;
	var y3_1_9;
	var z3_1_9;
	var x4_1_9;
	var y4_1_9;
	var z4_1_9;
	var x5_1_9;
	var y5_1_9;
	var z5_1_9;
	var x6_1_9;
	var y6_1_9;
	var z6_1_9;
	var x7_1_9;
	var y7_1_9;
	var z7_1_9;
	var x8_1_9;
	var y8_1_9;
	var z8_1_9;
	var x9_1_9;
	var y9_1_9;
	var z9_1_9;
	var x1_2_9;
	var y1_2_9;
	var z1_2_9;
	var x2_2_9;
	var y2_2_9;
	var z2_2_9;
	var x3_2_9;
	var y3_2_9;
	var z3_2_9;
	var x4_2_9;
	var y4_2_9;
	var z4_2_9;
	var x5_2_9;
	var y5_2_9;
	var z5_2_9;
	var x6_2_9;
	var y6_2_9;
	var z6_2_9;
	var x7_2_9;
	var y7_2_9;
	var z7_2_9;
	var x8_2_9;
	var y8_2_9;
	var z8_2_9;
	var x9_2_9;
	var y9_2_9;
	var z9_2_9;
	var x1_3_9;
	var y1_3_9;
	var z1_3_9;
	var x2_3_9;
	var y2_3_9;
	var z2_3_9;
	var x3_3_9;
	var y3_3_9;
	var z3_3_9;
	var x4_3_9;
	var y4_3_9;
	var z4_3_9;
	var x5_3_9;
	var y5_3_9;
	var z5_3_9;
	var x6_3_9;
	var y6_3_9;
	var z6_3_9;
	var x7_3_9;
	var y7_3_9;
	var z7_3_9;
	var x8_3_9;
	var y8_3_9;
	var z8_3_9;
	var x9_3_9;
	var y9_3_9;
	var z9_3_9;
	var x1_4_9;
	var y1_4_9;
	var z1_4_9;
	var x2_4_9;
	var y2_4_9;
	var z2_4_9;
	var x3_4_9;
	var y3_4_9;
	var z3_4_9;
	var x4_4_9;
	var y4_4_9;
	var z4_4_9;
	var x5_4_9;
	var y5_4_9;
	var z5_4_9;
	var x6_4_9;
	var y6_4_9;
	var z6_4_9;
	var x7_4_9;
	var y7_4_9;
	var z7_4_9;
	var x8_4_9;
	var y8_4_9;
	var z8_4_9;
	var x9_4_9;
	var y9_4_9;
	var z9_4_9;
	var x1_5_9;
	var y1_5_9;
	var z1_5_9;
	var x2_5_9;
	var y2_5_9;
	var z2_5_9;
	var x3_5_9;
	var y3_5_9;
	var z3_5_9;
	var x4_5_9;
	var y4_5_9;
	var z4_5_9;
	var x5_5_9;
	var y5_5_9;
	var z5_5_9;
	var x6_5_9;
	var y6_5_9;
	var z6_5_9;
	var x7_5_9;
	var y7_5_9;
	var z7_5_9;
	var x8_5_9;
	var y8_5_9;
	var z8_5_9;
	var x9_5_9;
	var y9_5_9;
	var z9_5_9;
	var x1_6_9;
	var y1_6_9;
	var z1_6_9;
	var x2_6_9;
	var y2_6_9;
	var z2_6_9;
	var x3_6_9;
	var y3_6_9;
	var z3_6_9;
	var x4_6_9;
	var y4_6_9;
	var z4_6_9;
	var x5_6_9;
	var y5_6_9;
	var z5_6_9;
	var x6_6_9;
	var y6_6_9;
	var z6_6_9;
	var x7_6_9;
	var y7_6_9;
	var z7_6_9;
	var x8_6_9;
	var y8_6_9;
	var z8_6_9;
	var x9_6_9;
	var y9_6_9;
	var z9_6_9;
	var x1_7_9;
	var y1_7_9;
	var z1_7_9;
	var x2_7_9;
	var y2_7_9;
	var z2_7_9;
	var x3_7_9;
	var y3_7_9;
	var z3_7_9;
	var x4_7_9;
	var y4_7_9;
	var z4_7_9;
	var x5_7_9;
	var y5_7_9;
	var z5_7_9;
	var x6_7_9;
	var y6_7_9;
	var z6_7_9;
	var x7_7_9;
	var y7_7_9;
	var z7_7_9;
	var x8_7_9;
	var y8_7_9;
	var z8_7_9;
	var x9_7_9;
	var y9_7_9;
	var z9_7_9;
	var x1_8_9;
	var y1_8_9;
	var z1_8_9;
	var x2_8_9;
	var y2_8_9;
	var z2_8_9;
	var x3_8_9;
	var y3_8_9;
	var z3_8_9;
	var x4_8_9;
	var y4_8_9;
	var z4_8_9;
	var x5_8_9;
	var y5_8_9;
	var z5_8_9;
	var x6_8_9;
	var y6_8_9;
	var z6_8_9;
	var x7_8_9;
	var y7_8_9;
	var z7_8_9;
	var x8_8_9;
	var y8_8_9;
	var z8_8_9;
	var x9_8_9;
	var y9_8_9;
	var z9_8_9;
	var x1_9_9;
	var y1_9_9;
	var z1_9_9;
	var x2_9_9;
	var y2_9_9;
	var z2_9_9;
	var x3_9_9;
	var y3_9_9;
	var z3_9_9;
	var x4_9_9;
	var y4_9_9;
	var z4_9_9;
	var x5_9_9;
	var y5_9_9;
	var z5_9_9;
	var x6_9_9;
	var y6_9_9;
	var z6_9_9;
	var x7_9_9;
	var y7_9_9;
	var z7_9_9;
	var x8_9_9;
	var y8_9_9;
	var z8_9_9;
	var x9_9_9;
	var y9_9_9;
	var z9_9_9;
	var y10_1_1;
	var z10_1_1;
	var y10_2_1;
	var z10_2_1;
	var y10_3_1;
	var z10_3_1;
	var y10_4_1;
	var z10_4_1;
	var y10_5_1;
	var z10_5_1;
	var y10_6_1;
	var z10_6_1;
	var y10_7_1;
	var z10_7_1;
	var y10_8_1;
	var z10_8_1;
	var y10_9_1;
	var z10_9_1;
	var y10_1_2;
	var z10_1_2;
	var y10_2_2;
	var z10_2_2;
	var y10_3_2;
	var z10_3_2;
	var y10_4_2;
	var z10_4_2;
	var y10_5_2;
	var z10_5_2;
	var y10_6_2;
	var z10_6_2;
	var y10_7_2;
	var z10_7_2;
	var y10_8_2;
	var z10_8_2;
	var y10_9_2;
	var z10_9_2;
	var y10_1_3;
	var z10_1_3;
	var y10_2_3;
	var z10_2_3;
	var y10_3_3;
	var z10_3_3;
	var y10_4_3;
	var z10_4_3;
	var y10_5_3;
	var z10_5_3;
	var y10_6_3;
	var z10_6_3;
	var y10_7_3;
	var z10_7_3;
	var y10_8_3;
	var z10_8_3;
	var y10_9_3;
	var z10_9_3;
	var y10_1_4;
	var z10_1_4;
	var y10_2_4;
	var z10_2_4;
	var y10_3_4;
	var z10_3_4;
	var y10_4_4;
	var z10_4_4;
	var y10_5_4;
	var z10_5_4;
	var y10_6_4;
	var z10_6_4;
	var y10_7_4;
	var z10_7_4;
	var y10_8_4;
	var z10_8_4;
	var y10_9_4;
	var z10_9_4;
	var y10_1_5;
	var z10_1_5;
	var y10_2_5;
	var z10_2_5;
	var y10_3_5;
	var z10_3_5;
	var y10_4_5;
	var z10_4_5;
	var y10_5_5;
	var z10_5_5;
	var y10_6_5;
	var z10_6_5;
	var y10_7_5;
	var z10_7_5;
	var y10_8_5;
	var z10_8_5;
	var y10_9_5;
	var z10_9_5;
	var y10_1_6;
	var z10_1_6;
	var y10_2_6;
	var z10_2_6;
	var y10_3_6;
	var z10_3_6;
	var y10_4_6;
	var z10_4_6;
	var y10_5_6;
	var z10_5_6;
	var y10_6_6;
	var z10_6_6;
	var y10_7_6;
	var z10_7_6;
	var y10_8_6;
	var z10_8_6;
	var y10_9_6;
	var z10_9_6;
	var y10_1_7;
	var z10_1_7;
	var y10_2_7;
	var z10_2_7;
	var y10_3_7;
	var z10_3_7;
	var y10_4_7;
	var z10_4_7;
	var y10_5_7;
	var z10_5_7;
	var y10_6_7;
	var z10_6_7;
	var y10_7_7;
	var z10_7_7;
	var y10_8_7;
	var z10_8_7;
	var y10_9_7;
	var z10_9_7;
	var y10_1_8;
	var z10_1_8;
	var y10_2_8;
	var z10_2_8;
	var y10_3_8;
	var z10_3_8;
	var y10_4_8;
	var z10_4_8;
	var y10_5_8;
	var z10_5_8;
	var y10_6_8;
	var z10_6_8;
	var y10_7_8;
	var z10_7_8;
	var y10_8_8;
	var z10_8_8;
	var y10_9_8;
	var z10_9_8;
	var y10_1_9;
	var z10_1_9;
	var y10_2_9;
	var z10_2_9;
	var y10_3_9;
	var z10_3_9;
	var y10_4_9;
	var z10_4_9;
	var y10_5_9;
	var z10_5_9;
	var y10_6_9;
	var z10_6_9;
	var y10_7_9;
	var z10_7_9;
	var y10_8_9;
	var z10_8_9;
	var y10_9_9;
	var z10_9_9;
	var x1_10_1;
	var z1_10_1;
	var x2_10_1;
	var z2_10_1;
	var x3_10_1;
	var z3_10_1;
	var x4_10_1;
	var z4_10_1;
	var x5_10_1;
	var z5_10_1;
	var x6_10_1;
	var z6_10_1;
	var x7_10_1;
	var z7_10_1;
	var x8_10_1;
	var z8_10_1;
	var x9_10_1;
	var z9_10_1;
	var x1_10_2;
	var z1_10_2;
	var x2_10_2;
	var z2_10_2;
	var x3_10_2;
	var z3_10_2;
	var x4_10_2;
	var z4_10_2;
	var x5_10_2;
	var z5_10_2;
	var x6_10_2;
	var z6_10_2;
	var x7_10_2;
	var z7_10_2;
	var x8_10_2;
	var z8_10_2;
	var x9_10_2;
	var z9_10_2;
	var x1_10_3;
	var z1_10_3;
	var x2_10_3;
	var z2_10_3;
	var x3_10_3;
	var z3_10_3;
	var x4_10_3;
	var z4_10_3;
	var x5_10_3;
	var z5_10_3;
	var x6_10_3;
	var z6_10_3;
	var x7_10_3;
	var z7_10_3;
	var x8_10_3;
	var z8_10_3;
	var x9_10_3;
	var z9_10_3;
	var x1_10_4;
	var z1_10_4;
	var x2_10_4;
	var z2_10_4;
	var x3_10_4;
	var z3_10_4;
	var x4_10_4;
	var z4_10_4;
	var x5_10_4;
	var z5_10_4;
	var x6_10_4;
	var z6_10_4;
	var x7_10_4;
	var z7_10_4;
	var x8_10_4;
	var z8_10_4;
	var x9_10_4;
	var z9_10_4;
	var x1_10_5;
	var z1_10_5;
	var x2_10_5;
	var z2_10_5;
	var x3_10_5;
	var z3_10_5;
	var x4_10_5;
	var z4_10_5;
	var x5_10_5;
	var z5_10_5;
	var x6_10_5;
	var z6_10_5;
	var x7_10_5;
	var z7_10_5;
	var x8_10_5;
	var z8_10_5;
	var x9_10_5;
	var z9_10_5;
	var x1_10_6;
	var z1_10_6;
	var x2_10_6;
	var z2_10_6;
	var x3_10_6;
	var z3_10_6;
	var x4_10_6;
	var z4_10_6;
	var x5_10_6;
	var z5_10_6;
	var x6_10_6;
	var z6_10_6;
	var x7_10_6;
	var z7_10_6;
	var x8_10_6;
	var z8_10_6;
	var x9_10_6;
	var z9_10_6;
	var x1_10_7;
	var z1_10_7;
	var x2_10_7;
	var z2_10_7;
	var x3_10_7;
	var z3_10_7;
	var x4_10_7;
	var z4_10_7;
	var x5_10_7;
	var z5_10_7;
	var x6_10_7;
	var z6_10_7;
	var x7_10_7;
	var z7_10_7;
	var x8_10_7;
	var z8_10_7;
	var x9_10_7;
	var z9_10_7;
	var x1_10_8;
	var z1_10_8;
	var x2_10_8;
	var z2_10_8;
	var x3_10_8;
	var z3_10_8;
	var x4_10_8;
	var z4_10_8;
	var x5_10_8;
	var z5_10_8;
	var x6_10_8;
	var z6_10_8;
	var x7_10_8;
	var z7_10_8;
	var x8_10_8;
	var z8_10_8;
	var x9_10_8;
	var z9_10_8;
	var x1_10_9;
	var z1_10_9;
	var x2_10_9;
	var z2_10_9;
	var x3_10_9;
	var z3_10_9;
	var x4_10_9;
	var z4_10_9;
	var x5_10_9;
	var z5_10_9;
	var x6_10_9;
	var z6_10_9;
	var x7_10_9;
	var z7_10_9;
	var x8_10_9;
	var z8_10_9;
	var x9_10_9;
	var z9_10_9;
	var x1_1_10;
	var y1_1_10;
	var x2_1_10;
	var y2_1_10;
	var x3_1_10;
	var y3_1_10;
	var x4_1_10;
	var y4_1_10;
	var x5_1_10;
	var y5_1_10;
	var x6_1_10;
	var y6_1_10;
	var x7_1_10;
	var y7_1_10;
	var x8_1_10;
	var y8_1_10;
	var x9_1_10;
	var y9_1_10;
	var x1_2_10;
	var y1_2_10;
	var x2_2_10;
	var y2_2_10;
	var x3_2_10;
	var y3_2_10;
	var x4_2_10;
	var y4_2_10;
	var x5_2_10;
	var y5_2_10;
	var x6_2_10;
	var y6_2_10;
	var x7_2_10;
	var y7_2_10;
	var x8_2_10;
	var y8_2_10;
	var x9_2_10;
	var y9_2_10;
	var x1_3_10;
	var y1_3_10;
	var x2_3_10;
	var y2_3_10;
	var x3_3_10;
	var y3_3_10;
	var x4_3_10;
	var y4_3_10;
	var x5_3_10;
	var y5_3_10;
	var x6_3_10;
	var y6_3_10;
	var x7_3_10;
	var y7_3_10;
	var x8_3_10;
	var y8_3_10;
	var x9_3_10;
	var y9_3_10;
	var x1_4_10;
	var y1_4_10;
	var x2_4_10;
	var y2_4_10;
	var x3_4_10;
	var y3_4_10;
	var x4_4_10;
	var y4_4_10;
	var x5_4_10;
	var y5_4_10;
	var x6_4_10;
	var y6_4_10;
	var x7_4_10;
	var y7_4_10;
	var x8_4_10;
	var y8_4_10;
	var x9_4_10;
	var y9_4_10;
	var x1_5_10;
	var y1_5_10;
	var x2_5_10;
	var y2_5_10;
	var x3_5_10;
	var y3_5_10;
	var x4_5_10;
	var y4_5_10;
	var x5_5_10;
	var y5_5_10;
	var x6_5_10;
	var y6_5_10;
	var x7_5_10;
	var y7_5_10;
	var x8_5_10;
	var y8_5_10;
	var x9_5_10;
	var y9_5_10;
	var x1_6_10;
	var y1_6_10;
	var x2_6_10;
	var y2_6_10;
	var x3_6_10;
	var y3_6_10;
	var x4_6_10;
	var y4_6_10;
	var x5_6_10;
	var y5_6_10;
	var x6_6_10;
	var y6_6_10;
	var x7_6_10;
	var y7_6_10;
	var x8_6_10;
	var y8_6_10;
	var x9_6_10;
	var y9_6_10;
	var x1_7_10;
	var y1_7_10;
	var x2_7_10;
	var y2_7_10;
	var x3_7_10;
	var y3_7_10;
	var x4_7_10;
	var y4_7_10;
	var x5_7_10;
	var y5_7_10;
	var x6_7_10;
	var y6_7_10;
	var x7_7_10;
	var y7_7_10;
	var x8_7_10;
	var y8_7_10;
	var x9_7_10;
	var y9_7_10;
	var x1_8_10;
	var y1_8_10;
	var x2_8_10;
	var y2_8_10;
	var x3_8_10;
	var y3_8_10;
	var x4_8_10;
	var y4_8_10;
	var x5_8_10;
	var y5_8_10;
	var x6_8_10;
	var y6_8_10;
	var x7_8_10;
	var y7_8_10;
	var x8_8_10;
	var y8_8_10;
	var x9_8_10;
	var y9_8_10;
	var x1_9_10;
	var y1_9_10;
	var x2_9_10;
	var y2_9_10;
	var x3_9_10;
	var y3_9_10;
	var x4_9_10;
	var y4_9_10;
	var x5_9_10;
	var y5_9_10;
	var x6_9_10;
	var y6_9_10;
	var x7_9_10;
	var y7_9_10;
	var x8_9_10;
	var y8_9_10;
	var x9_9_10;
	var y9_9_10;
	var y0_1_1;
	var z0_1_1;
	var y11_1_1;
	var z11_1_1;
	var y0_2_1;
	var z0_2_1;
	var y11_2_1;
	var z11_2_1;
	var y0_3_1;
	var z0_3_1;
	var y11_3_1;
	var z11_3_1;
	var y0_4_1;
	var z0_4_1;
	var y11_4_1;
	var z11_4_1;
	var y0_5_1;
	var z0_5_1;
	var y11_5_1;
	var z11_5_1;
	var y0_6_1;
	var z0_6_1;
	var y11_6_1;
	var z11_6_1;
	var y0_7_1;
	var z0_7_1;
	var y11_7_1;
	var z11_7_1;
	var y0_8_1;
	var z0_8_1;
	var y11_8_1;
	var z11_8_1;
	var y0_9_1;
	var z0_9_1;
	var y11_9_1;
	var z11_9_1;
	var y0_10_1;
	var z0_10_1;
	var y11_10_1;
	var z11_10_1;
	var y0_1_2;
	var z0_1_2;
	var y11_1_2;
	var z11_1_2;
	var y0_2_2;
	var z0_2_2;
	var y11_2_2;
	var z11_2_2;
	var y0_3_2;
	var z0_3_2;
	var y11_3_2;
	var z11_3_2;
	var y0_4_2;
	var z0_4_2;
	var y11_4_2;
	var z11_4_2;
	var y0_5_2;
	var z0_5_2;
	var y11_5_2;
	var z11_5_2;
	var y0_6_2;
	var z0_6_2;
	var y11_6_2;
	var z11_6_2;
	var y0_7_2;
	var z0_7_2;
	var y11_7_2;
	var z11_7_2;
	var y0_8_2;
	var z0_8_2;
	var y11_8_2;
	var z11_8_2;
	var y0_9_2;
	var z0_9_2;
	var y11_9_2;
	var z11_9_2;
	var y0_10_2;
	var z0_10_2;
	var y11_10_2;
	var z11_10_2;
	var y0_1_3;
	var z0_1_3;
	var y11_1_3;
	var z11_1_3;
	var y0_2_3;
	var z0_2_3;
	var y11_2_3;
	var z11_2_3;
	var y0_3_3;
	var z0_3_3;
	var y11_3_3;
	var z11_3_3;
	var y0_4_3;
	var z0_4_3;
	var y11_4_3;
	var z11_4_3;
	var y0_5_3;
	var z0_5_3;
	var y11_5_3;
	var z11_5_3;
	var y0_6_3;
	var z0_6_3;
	var y11_6_3;
	var z11_6_3;
	var y0_7_3;
	var z0_7_3;
	var y11_7_3;
	var z11_7_3;
	var y0_8_3;
	var z0_8_3;
	var y11_8_3;
	var z11_8_3;
	var y0_9_3;
	var z0_9_3;
	var y11_9_3;
	var z11_9_3;
	var y0_10_3;
	var z0_10_3;
	var y11_10_3;
	var z11_10_3;
	var y0_1_4;
	var z0_1_4;
	var y11_1_4;
	var z11_1_4;
	var y0_2_4;
	var z0_2_4;
	var y11_2_4;
	var z11_2_4;
	var y0_3_4;
	var z0_3_4;
	var y11_3_4;
	var z11_3_4;
	var y0_4_4;
	var z0_4_4;
	var y11_4_4;
	var z11_4_4;
	var y0_5_4;
	var z0_5_4;
	var y11_5_4;
	var z11_5_4;
	var y0_6_4;
	var z0_6_4;
	var y11_6_4;
	var z11_6_4;
	var y0_7_4;
	var z0_7_4;
	var y11_7_4;
	var z11_7_4;
	var y0_8_4;
	var z0_8_4;
	var y11_8_4;
	var z11_8_4;
	var y0_9_4;
	var z0_9_4;
	var y11_9_4;
	var z11_9_4;
	var y0_10_4;
	var z0_10_4;
	var y11_10_4;
	var z11_10_4;
	var y0_1_5;
	var z0_1_5;
	var y11_1_5;
	var z11_1_5;
	var y0_2_5;
	var z0_2_5;
	var y11_2_5;
	var z11_2_5;
	var y0_3_5;
	var z0_3_5;
	var y11_3_5;
	var z11_3_5;
	var y0_4_5;
	var z0_4_5;
	var y11_4_5;
	var z11_4_5;
	var y0_5_5;
	var z0_5_5;
	var y11_5_5;
	var z11_5_5;
	var y0_6_5;
	var z0_6_5;
	var y11_6_5;
	var z11_6_5;
	var y0_7_5;
	var z0_7_5;
	var y11_7_5;
	var z11_7_5;
	var y0_8_5;
	var z0_8_5;
	var y11_8_5;
	var z11_8_5;
	var y0_9_5;
	var z0_9_5;
	var y11_9_5;
	var z11_9_5;
	var y0_10_5;
	var z0_10_5;
	var y11_10_5;
	var z11_10_5;
	var y0_1_6;
	var z0_1_6;
	var y11_1_6;
	var z11_1_6;
	var y0_2_6;
	var z0_2_6;
	var y11_2_6;
	var z11_2_6;
	var y0_3_6;
	var z0_3_6;
	var y11_3_6;
	var z11_3_6;
	var y0_4_6;
	var z0_4_6;
	var y11_4_6;
	var z11_4_6;
	var y0_5_6;
	var z0_5_6;
	var y11_5_6;
	var z11_5_6;
	var y0_6_6;
	var z0_6_6;
	var y11_6_6;
	var z11_6_6;
	var y0_7_6;
	var z0_7_6;
	var y11_7_6;
	var z11_7_6;
	var y0_8_6;
	var z0_8_6;
	var y11_8_6;
	var z11_8_6;
	var y0_9_6;
	var z0_9_6;
	var y11_9_6;
	var z11_9_6;
	var y0_10_6;
	var z0_10_6;
	var y11_10_6;
	var z11_10_6;
	var y0_1_7;
	var z0_1_7;
	var y11_1_7;
	var z11_1_7;
	var y0_2_7;
	var z0_2_7;
	var y11_2_7;
	var z11_2_7;
	var y0_3_7;
	var z0_3_7;
	var y11_3_7;
	var z11_3_7;
	var y0_4_7;
	var z0_4_7;
	var y11_4_7;
	var z11_4_7;
	var y0_5_7;
	var z0_5_7;
	var y11_5_7;
	var z11_5_7;
	var y0_6_7;
	var z0_6_7;
	var y11_6_7;
	var z11_6_7;
	var y0_7_7;
	var z0_7_7;
	var y11_7_7;
	var z11_7_7;
	var y0_8_7;
	var z0_8_7;
	var y11_8_7;
	var z11_8_7;
	var y0_9_7;
	var z0_9_7;
	var y11_9_7;
	var z11_9_7;
	var y0_10_7;
	var z0_10_7;
	var y11_10_7;
	var z11_10_7;
	var y0_1_8;
	var z0_1_8;
	var y11_1_8;
	var z11_1_8;
	var y0_2_8;
	var z0_2_8;
	var y11_2_8;
	var z11_2_8;
	var y0_3_8;
	var z0_3_8;
	var y11_3_8;
	var z11_3_8;
	var y0_4_8;
	var z0_4_8;
	var y11_4_8;
	var z11_4_8;
	var y0_5_8;
	var z0_5_8;
	var y11_5_8;
	var z11_5_8;
	var y0_6_8;
	var z0_6_8;
	var y11_6_8;
	var z11_6_8;
	var y0_7_8;
	var z0_7_8;
	var y11_7_8;
	var z11_7_8;
	var y0_8_8;
	var z0_8_8;
	var y11_8_8;
	var z11_8_8;
	var y0_9_8;
	var z0_9_8;
	var y11_9_8;
	var z11_9_8;
	var y0_10_8;
	var z0_10_8;
	var y11_10_8;
	var z11_10_8;
	var y0_1_9;
	var z0_1_9;
	var y11_1_9;
	var z11_1_9;
	var y0_2_9;
	var z0_2_9;
	var y11_2_9;
	var z11_2_9;
	var y0_3_9;
	var z0_3_9;
	var y11_3_9;
	var z11_3_9;
	var y0_4_9;
	var z0_4_9;
	var y11_4_9;
	var z11_4_9;
	var y0_5_9;
	var z0_5_9;
	var y11_5_9;
	var z11_5_9;
	var y0_6_9;
	var z0_6_9;
	var y11_6_9;
	var z11_6_9;
	var y0_7_9;
	var z0_7_9;
	var y11_7_9;
	var z11_7_9;
	var y0_8_9;
	var z0_8_9;
	var y11_8_9;
	var z11_8_9;
	var y0_9_9;
	var z0_9_9;
	var y11_9_9;
	var z11_9_9;
	var y0_10_9;
	var z0_10_9;
	var y11_10_9;
	var z11_10_9;
	var y0_1_10;
	var z0_1_10;
	var y11_1_10;
	var z11_1_10;
	var y0_2_10;
	var z0_2_10;
	var y11_2_10;
	var z11_2_10;
	var y0_3_10;
	var z0_3_10;
	var y11_3_10;
	var z11_3_10;
	var y0_4_10;
	var z0_4_10;
	var y11_4_10;
	var z11_4_10;
	var y0_5_10;
	var z0_5_10;
	var y11_5_10;
	var z11_5_10;
	var y0_6_10;
	var z0_6_10;
	var y11_6_10;
	var z11_6_10;
	var y0_7_10;
	var z0_7_10;
	var y11_7_10;
	var z11_7_10;
	var y0_8_10;
	var z0_8_10;
	var y11_8_10;
	var z11_8_10;
	var y0_9_10;
	var z0_9_10;
	var y11_9_10;
	var z11_9_10;
	var y0_10_10;
	var z0_10_10;
	var y11_10_10;
	var z11_10_10;
	var x1_0_1;
	var z1_0_1;
	var x1_11_1;
	var z1_11_1;
	var x2_0_1;
	var z2_0_1;
	var x2_11_1;
	var z2_11_1;
	var x3_0_1;
	var z3_0_1;
	var x3_11_1;
	var z3_11_1;
	var x4_0_1;
	var z4_0_1;
	var x4_11_1;
	var z4_11_1;
	var x5_0_1;
	var z5_0_1;
	var x5_11_1;
	var z5_11_1;
	var x6_0_1;
	var z6_0_1;
	var x6_11_1;
	var z6_11_1;
	var x7_0_1;
	var z7_0_1;
	var x7_11_1;
	var z7_11_1;
	var x8_0_1;
	var z8_0_1;
	var x8_11_1;
	var z8_11_1;
	var x9_0_1;
	var z9_0_1;
	var x9_11_1;
	var z9_11_1;
	var x10_0_1;
	var z10_0_1;
	var x10_11_1;
	var z10_11_1;
	var x1_0_2;
	var z1_0_2;
	var x1_11_2;
	var z1_11_2;
	var x2_0_2;
	var z2_0_2;
	var x2_11_2;
	var z2_11_2;
	var x3_0_2;
	var z3_0_2;
	var x3_11_2;
	var z3_11_2;
	var x4_0_2;
	var z4_0_2;
	var x4_11_2;
	var z4_11_2;
	var x5_0_2;
	var z5_0_2;
	var x5_11_2;
	var z5_11_2;
	var x6_0_2;
	var z6_0_2;
	var x6_11_2;
	var z6_11_2;
	var x7_0_2;
	var z7_0_2;
	var x7_11_2;
	var z7_11_2;
	var x8_0_2;
	var z8_0_2;
	var x8_11_2;
	var z8_11_2;
	var x9_0_2;
	var z9_0_2;
	var x9_11_2;
	var z9_11_2;
	var x10_0_2;
	var z10_0_2;
	var x10_11_2;
	var z10_11_2;
	var x1_0_3;
	var z1_0_3;
	var x1_11_3;
	var z1_11_3;
	var x2_0_3;
	var z2_0_3;
	var x2_11_3;
	var z2_11_3;
	var x3_0_3;
	var z3_0_3;
	var x3_11_3;
	var z3_11_3;
	var x4_0_3;
	var z4_0_3;
	var x4_11_3;
	var z4_11_3;
	var x5_0_3;
	var z5_0_3;
	var x5_11_3;
	var z5_11_3;
	var x6_0_3;
	var z6_0_3;
	var x6_11_3;
	var z6_11_3;
	var x7_0_3;
	var z7_0_3;
	var x7_11_3;
	var z7_11_3;
	var x8_0_3;
	var z8_0_3;
	var x8_11_3;
	var z8_11_3;
	var x9_0_3;
	var z9_0_3;
	var x9_11_3;
	var z9_11_3;
	var x10_0_3;
	var z10_0_3;
	var x10_11_3;
	var z10_11_3;
	var x1_0_4;
	var z1_0_4;
	var x1_11_4;
	var z1_11_4;
	var x2_0_4;
	var z2_0_4;
	var x2_11_4;
	var z2_11_4;
	var x3_0_4;
	var z3_0_4;
	var x3_11_4;
	var z3_11_4;
	var x4_0_4;
	var z4_0_4;
	var x4_11_4;
	var z4_11_4;
	var x5_0_4;
	var z5_0_4;
	var x5_11_4;
	var z5_11_4;
	var x6_0_4;
	var z6_0_4;
	var x6_11_4;
	var z6_11_4;
	var x7_0_4;
	var z7_0_4;
	var x7_11_4;
	var z7_11_4;
	var x8_0_4;
	var z8_0_4;
	var x8_11_4;
	var z8_11_4;
	var x9_0_4;
	var z9_0_4;
	var x9_11_4;
	var z9_11_4;
	var x10_0_4;
	var z10_0_4;
	var x10_11_4;
	var z10_11_4;
	var x1_0_5;
	var z1_0_5;
	var x1_11_5;
	var z1_11_5;
	var x2_0_5;
	var z2_0_5;
	var x2_11_5;
	var z2_11_5;
	var x3_0_5;
	var z3_0_5;
	var x3_11_5;
	var z3_11_5;
	var x4_0_5;
	var z4_0_5;
	var x4_11_5;
	var z4_11_5;
	var x5_0_5;
	var z5_0_5;
	var x5_11_5;
	var z5_11_5;
	var x6_0_5;
	var z6_0_5;
	var x6_11_5;
	var z6_11_5;
	var x7_0_5;
	var z7_0_5;
	var x7_11_5;
	var z7_11_5;
	var x8_0_5;
	var z8_0_5;
	var x8_11_5;
	var z8_11_5;
	var x9_0_5;
	var z9_0_5;
	var x9_11_5;
	var z9_11_5;
	var x10_0_5;
	var z10_0_5;
	var x10_11_5;
	var z10_11_5;
	var x1_0_6;
	var z1_0_6;
	var x1_11_6;
	var z1_11_6;
	var x2_0_6;
	var z2_0_6;
	var x2_11_6;
	var z2_11_6;
	var x3_0_6;
	var z3_0_6;
	var x3_11_6;
	var z3_11_6;
	var x4_0_6;
	var z4_0_6;
	var x4_11_6;
	var z4_11_6;
	var x5_0_6;
	var z5_0_6;
	var x5_11_6;
	var z5_11_6;
	var x6_0_6;
	var z6_0_6;
	var x6_11_6;
	var z6_11_6;
	var x7_0_6;
	var z7_0_6;
	var x7_11_6;
	var z7_11_6;
	var x8_0_6;
	var z8_0_6;
	var x8_11_6;
	var z8_11_6;
	var x9_0_6;
	var z9_0_6;
	var x9_11_6;
	var z9_11_6;
	var x10_0_6;
	var z10_0_6;
	var x10_11_6;
	var z10_11_6;
	var x1_0_7;
	var z1_0_7;
	var x1_11_7;
	var z1_11_7;
	var x2_0_7;
	var z2_0_7;
	var x2_11_7;
	var z2_11_7;
	var x3_0_7;
	var z3_0_7;
	var x3_11_7;
	var z3_11_7;
	var x4_0_7;
	var z4_0_7;
	var x4_11_7;
	var z4_11_7;
	var x5_0_7;
	var z5_0_7;
	var x5_11_7;
	var z5_11_7;
	var x6_0_7;
	var z6_0_7;
	var x6_11_7;
	var z6_11_7;
	var x7_0_7;
	var z7_0_7;
	var x7_11_7;
	var z7_11_7;
	var x8_0_7;
	var z8_0_7;
	var x8_11_7;
	var z8_11_7;
	var x9_0_7;
	var z9_0_7;
	var x9_11_7;
	var z9_11_7;
	var x10_0_7;
	var z10_0_7;
	var x10_11_7;
	var z10_11_7;
	var x1_0_8;
	var z1_0_8;
	var x1_11_8;
	var z1_11_8;
	var x2_0_8;
	var z2_0_8;
	var x2_11_8;
	var z2_11_8;
	var x3_0_8;
	var z3_0_8;
	var x3_11_8;
	var z3_11_8;
	var x4_0_8;
	var z4_0_8;
	var x4_11_8;
	var z4_11_8;
	var x5_0_8;
	var z5_0_8;
	var x5_11_8;
	var z5_11_8;
	var x6_0_8;
	var z6_0_8;
	var x6_11_8;
	var z6_11_8;
	var x7_0_8;
	var z7_0_8;
	var x7_11_8;
	var z7_11_8;
	var x8_0_8;
	var z8_0_8;
	var x8_11_8;
	var z8_11_8;
	var x9_0_8;
	var z9_0_8;
	var x9_11_8;
	var z9_11_8;
	var x10_0_8;
	var z10_0_8;
	var x10_11_8;
	var z10_11_8;
	var x1_0_9;
	var z1_0_9;
	var x1_11_9;
	var z1_11_9;
	var x2_0_9;
	var z2_0_9;
	var x2_11_9;
	var z2_11_9;
	var x3_0_9;
	var z3_0_9;
	var x3_11_9;
	var z3_11_9;
	var x4_0_9;
	var z4_0_9;
	var x4_11_9;
	var z4_11_9;
	var x5_0_9;
	var z5_0_9;
	var x5_11_9;
	var z5_11_9;
	var x6_0_9;
	var z6_0_9;
	var x6_11_9;
	var z6_11_9;
	var x7_0_9;
	var z7_0_9;
	var x7_11_9;
	var z7_11_9;
	var x8_0_9;
	var z8_0_9;
	var x8_11_9;
	var z8_11_9;
	var x9_0_9;
	var z9_0_9;
	var x9_11_9;
	var z9_11_9;
	var x10_0_9;
	var z10_0_9;
	var x10_11_9;
	var z10_11_9;
	var x1_0_10;
	var z1_0_10;
	var x1_11_10;
	var z1_11_10;
	var x2_0_10;
	var z2_0_10;
	var x2_11_10;
	var z2_11_10;
	var x3_0_10;
	var z3_0_10;
	var x3_11_10;
	var z3_11_10;
	var x4_0_10;
	var z4_0_10;
	var x4_11_10;
	var z4_11_10;
	var x5_0_10;
	var z5_0_10;
	var x5_11_10;
	var z5_11_10;
	var x6_0_10;
	var z6_0_10;
	var x6_11_10;
	var z6_11_10;
	var x7_0_10;
	var z7_0_10;
	var x7_11_10;
	var z7_11_10;
	var x8_0_10;
	var z8_0_10;
	var x8_11_10;
	var z8_11_10;
	var x9_0_10;
	var z9_0_10;
	var x9_11_10;
	var z9_11_10;
	var x10_0_10;
	var z10_0_10;
	var x10_11_10;
	var z10_11_10;
	var x1_1_0;
	var y1_1_0;
	var x1_1_11;
	var y1_1_11;
	var x2_1_0;
	var y2_1_0;
	var x2_1_11;
	var y2_1_11;
	var x3_1_0;
	var y3_1_0;
	var x3_1_11;
	var y3_1_11;
	var x4_1_0;
	var y4_1_0;
	var x4_1_11;
	var y4_1_11;
	var x5_1_0;
	var y5_1_0;
	var x5_1_11;
	var y5_1_11;
	var x6_1_0;
	var y6_1_0;
	var x6_1_11;
	var y6_1_11;
	var x7_1_0;
	var y7_1_0;
	var x7_1_11;
	var y7_1_11;
	var x8_1_0;
	var y8_1_0;
	var x8_1_11;
	var y8_1_11;
	var x9_1_0;
	var y9_1_0;
	var x9_1_11;
	var y9_1_11;
	var x10_1_0;
	var y10_1_0;
	var x10_1_11;
	var y10_1_11;
	var x1_2_0;
	var y1_2_0;
	var x1_2_11;
	var y1_2_11;
	var x2_2_0;
	var y2_2_0;
	var x2_2_11;
	var y2_2_11;
	var x3_2_0;
	var y3_2_0;
	var x3_2_11;
	var y3_2_11;
	var x4_2_0;
	var y4_2_0;
	var x4_2_11;
	var y4_2_11;
	var x5_2_0;
	var y5_2_0;
	var x5_2_11;
	var y5_2_11;
	var x6_2_0;
	var y6_2_0;
	var x6_2_11;
	var y6_2_11;
	var x7_2_0;
	var y7_2_0;
	var x7_2_11;
	var y7_2_11;
	var x8_2_0;
	var y8_2_0;
	var x8_2_11;
	var y8_2_11;
	var x9_2_0;
	var y9_2_0;
	var x9_2_11;
	var y9_2_11;
	var x10_2_0;
	var y10_2_0;
	var x10_2_11;
	var y10_2_11;
	var x1_3_0;
	var y1_3_0;
	var x1_3_11;
	var y1_3_11;
	var x2_3_0;
	var y2_3_0;
	var x2_3_11;
	var y2_3_11;
	var x3_3_0;
	var y3_3_0;
	var x3_3_11;
	var y3_3_11;
	var x4_3_0;
	var y4_3_0;
	var x4_3_11;
	var y4_3_11;
	var x5_3_0;
	var y5_3_0;
	var x5_3_11;
	var y5_3_11;
	var x6_3_0;
	var y6_3_0;
	var x6_3_11;
	var y6_3_11;
	var x7_3_0;
	var y7_3_0;
	var x7_3_11;
	var y7_3_11;
	var x8_3_0;
	var y8_3_0;
	var x8_3_11;
	var y8_3_11;
	var x9_3_0;
	var y9_3_0;
	var x9_3_11;
	var y9_3_11;
	var x10_3_0;
	var y10_3_0;
	var x10_3_11;
	var y10_3_11;
	var x1_4_0;
	var y1_4_0;
	var x1_4_11;
	var y1_4_11;
	var x2_4_0;
	var y2_4_0;
	var x2_4_11;
	var y2_4_11;
	var x3_4_0;
	var y3_4_0;
	var x3_4_11;
	var y3_4_11;
	var x4_4_0;
	var y4_4_0;
	var x4_4_11;
	var y4_4_11;
	var x5_4_0;
	var y5_4_0;
	var x5_4_11;
	var y5_4_11;
	var x6_4_0;
	var y6_4_0;
	var x6_4_11;
	var y6_4_11;
	var x7_4_0;
	var y7_4_0;
	var x7_4_11;
	var y7_4_11;
	var x8_4_0;
	var y8_4_0;
	var x8_4_11;
	var y8_4_11;
	var x9_4_0;
	var y9_4_0;
	var x9_4_11;
	var y9_4_11;
	var x10_4_0;
	var y10_4_0;
	var x10_4_11;
	var y10_4_11;
	var x1_5_0;
	var y1_5_0;
	var x1_5_11;
	var y1_5_11;
	var x2_5_0;
	var y2_5_0;
	var x2_5_11;
	var y2_5_11;
	var x3_5_0;
	var y3_5_0;
	var x3_5_11;
	var y3_5_11;
	var x4_5_0;
	var y4_5_0;
	var x4_5_11;
	var y4_5_11;
	var x5_5_0;
	var y5_5_0;
	var x5_5_11;
	var y5_5_11;
	var x6_5_0;
	var y6_5_0;
	var x6_5_11;
	var y6_5_11;
	var x7_5_0;
	var y7_5_0;
	var x7_5_11;
	var y7_5_11;
	var x8_5_0;
	var y8_5_0;
	var x8_5_11;
	var y8_5_11;
	var x9_5_0;
	var y9_5_0;
	var x9_5_11;
	var y9_5_11;
	var x10_5_0;
	var y10_5_0;
	var x10_5_11;
	var y10_5_11;
	var x1_6_0;
	var y1_6_0;
	var x1_6_11;
	var y1_6_11;
	var x2_6_0;
	var y2_6_0;
	var x2_6_11;
	var y2_6_11;
	var x3_6_0;
	var y3_6_0;
	var x3_6_11;
	var y3_6_11;
	var x4_6_0;
	var y4_6_0;
	var x4_6_11;
	var y4_6_11;
	var x5_6_0;
	var y5_6_0;
	var x5_6_11;
	var y5_6_11;
	var x6_6_0;
	var y6_6_0;
	var x6_6_11;
	var y6_6_11;
	var x7_6_0;
	var y7_6_0;
	var x7_6_11;
	var y7_6_11;
	var x8_6_0;
	var y8_6_0;
	var x8_6_11;
	var y8_6_11;
	var x9_6_0;
	var y9_6_0;
	var x9_6_11;
	var y9_6_11;
	var x10_6_0;
	var y10_6_0;
	var x10_6_11;
	var y10_6_11;
	var x1_7_0;
	var y1_7_0;
	var x1_7_11;
	var y1_7_11;
	var x2_7_0;
	var y2_7_0;
	var x2_7_11;
	var y2_7_11;
	var x3_7_0;
	var y3_7_0;
	var x3_7_11;
	var y3_7_11;
	var x4_7_0;
	var y4_7_0;
	var x4_7_11;
	var y4_7_11;
	var x5_7_0;
	var y5_7_0;
	var x5_7_11;
	var y5_7_11;
	var x6_7_0;
	var y6_7_0;
	var x6_7_11;
	var y6_7_11;
	var x7_7_0;
	var y7_7_0;
	var x7_7_11;
	var y7_7_11;
	var x8_7_0;
	var y8_7_0;
	var x8_7_11;
	var y8_7_11;
	var x9_7_0;
	var y9_7_0;
	var x9_7_11;
	var y9_7_11;
	var x10_7_0;
	var y10_7_0;
	var x10_7_11;
	var y10_7_11;
	var x1_8_0;
	var y1_8_0;
	var x1_8_11;
	var y1_8_11;
	var x2_8_0;
	var y2_8_0;
	var x2_8_11;
	var y2_8_11;
	var x3_8_0;
	var y3_8_0;
	var x3_8_11;
	var y3_8_11;
	var x4_8_0;
	var y4_8_0;
	var x4_8_11;
	var y4_8_11;
	var x5_8_0;
	var y5_8_0;
	var x5_8_11;
	var y5_8_11;
	var x6_8_0;
	var y6_8_0;
	var x6_8_11;
	var y6_8_11;
	var x7_8_0;
	var y7_8_0;
	var x7_8_11;
	var y7_8_11;
	var x8_8_0;
	var y8_8_0;
	var x8_8_11;
	var y8_8_11;
	var x9_8_0;
	var y9_8_0;
	var x9_8_11;
	var y9_8_11;
	var x10_8_0;
	var y10_8_0;
	var x10_8_11;
	var y10_8_11;
	var x1_9_0;
	var y1_9_0;
	var x1_9_11;
	var y1_9_11;
	var x2_9_0;
	var y2_9_0;
	var x2_9_11;
	var y2_9_11;
	var x3_9_0;
	var y3_9_0;
	var x3_9_11;
	var y3_9_11;
	var x4_9_0;
	var y4_9_0;
	var x4_9_11;
	var y4_9_11;
	var x5_9_0;
	var y5_9_0;
	var x5_9_11;
	var y5_9_11;
	var x6_9_0;
	var y6_9_0;
	var x6_9_11;
	var y6_9_11;
	var x7_9_0;
	var y7_9_0;
	var x7_9_11;
	var y7_9_11;
	var x8_9_0;
	var y8_9_0;
	var x8_9_11;
	var y8_9_11;
	var x9_9_0;
	var y9_9_0;
	var x9_9_11;
	var y9_9_11;
	var x10_9_0;
	var y10_9_0;
	var x10_9_11;
	var y10_9_11;
	var x1_10_0;
	var y1_10_0;
	var x1_10_11;
	var y1_10_11;
	var x2_10_0;
	var y2_10_0;
	var x2_10_11;
	var y2_10_11;
	var x3_10_0;
	var y3_10_0;
	var x3_10_11;
	var y3_10_11;
	var x4_10_0;
	var y4_10_0;
	var x4_10_11;
	var y4_10_11;
	var x5_10_0;
	var y5_10_0;
	var x5_10_11;
	var y5_10_11;
	var x6_10_0;
	var y6_10_0;
	var x6_10_11;
	var y6_10_11;
	var x7_10_0;
	var y7_10_0;
	var x7_10_11;
	var y7_10_11;
	var x8_10_0;
	var y8_10_0;
	var x8_10_11;
	var y8_10_11;
	var x9_10_0;
	var y9_10_0;
	var x9_10_11;
	var y9_10_11;
	var x10_10_0;
	var y10_10_0;
	var x10_10_11;
	var y10_10_11;

minimize obj:
	5.0d-1*(x1_1_1 - 1.0)*(x1_1_1 - 1.0) + 5.0d-1*(y1_1_1 - 1.0)*(y1_1_1 - 1.0) + 
	5.0d-1*(z1_1_1 - 1.0)*(z1_1_1 - 1.0) + 5.0d-1*(x2_1_1 - 1.0)*(x2_1_1 - 1.0) + 
	5.0d-1*(y2_1_1 - 1.0)*(y2_1_1 - 1.0) + 5.0d-1*(z2_1_1 - 1.0)*(z2_1_1 - 1.0) + 
	5.0d-1*(x3_1_1 - 1.0)*(x3_1_1 - 1.0) + 5.0d-1*(y3_1_1 - 1.0)*(y3_1_1 - 1.0) + 
	5.0d-1*(z3_1_1 - 1.0)*(z3_1_1 - 1.0) + 5.0d-1*(x4_1_1 - 1.0)*(x4_1_1 - 1.0) + 
	5.0d-1*(y4_1_1 - 1.0)*(y4_1_1 - 1.0) + 5.0d-1*(z4_1_1 - 1.0)*(z4_1_1 - 1.0) + 
	5.0d-1*(x5_1_1 - 1.0)*(x5_1_1 - 1.0) + 5.0d-1*(y5_1_1 - 1.0)*(y5_1_1 - 1.0) + 
	5.0d-1*(z5_1_1 - 1.0)*(z5_1_1 - 1.0) + 5.0d-1*(x6_1_1 - 1.0)*(x6_1_1 - 1.0) + 
	5.0d-1*(y6_1_1 - 1.0)*(y6_1_1 - 1.0) + 5.0d-1*(z6_1_1 - 1.0)*(z6_1_1 - 1.0) + 
	5.0d-1*(x7_1_1 - 1.0)*(x7_1_1 - 1.0) + 5.0d-1*(y7_1_1 - 1.0)*(y7_1_1 - 1.0) + 
	5.0d-1*(z7_1_1 - 1.0)*(z7_1_1 - 1.0) + 5.0d-1*(x8_1_1 - 1.0)*(x8_1_1 - 1.0) + 
	5.0d-1*(y8_1_1 - 1.0)*(y8_1_1 - 1.0) + 5.0d-1*(z8_1_1 - 1.0)*(z8_1_1 - 1.0) + 
	5.0d-1*(x9_1_1 - 1.0)*(x9_1_1 - 1.0) + 5.0d-1*(y9_1_1 - 1.0)*(y9_1_1 - 1.0) + 
	5.0d-1*(z9_1_1 - 1.0)*(z9_1_1 - 1.0) + 5.0d-1*(x1_2_1 - 1.0)*(x1_2_1 - 1.0) + 
	5.0d-1*(y1_2_1 - 1.0)*(y1_2_1 - 1.0) + 5.0d-1*(z1_2_1 - 1.0)*(z1_2_1 - 1.0) + 
	5.0d-1*(x2_2_1 - 1.0)*(x2_2_1 - 1.0) + 5.0d-1*(y2_2_1 - 1.0)*(y2_2_1 - 1.0) + 
	5.0d-1*(z2_2_1 - 1.0)*(z2_2_1 - 1.0) + 5.0d-1*(x3_2_1 - 1.0)*(x3_2_1 - 1.0) + 
	5.0d-1*(y3_2_1 - 1.0)*(y3_2_1 - 1.0) + 5.0d-1*(z3_2_1 - 1.0)*(z3_2_1 - 1.0) + 
	5.0d-1*(x4_2_1 - 1.0)*(x4_2_1 - 1.0) + 5.0d-1*(y4_2_1 - 1.0)*(y4_2_1 - 1.0) + 
	5.0d-1*(z4_2_1 - 1.0)*(z4_2_1 - 1.0) + 5.0d-1*(x5_2_1 - 1.0)*(x5_2_1 - 1.0) + 
	5.0d-1*(y5_2_1 - 1.0)*(y5_2_1 - 1.0) + 5.0d-1*(z5_2_1 - 1.0)*(z5_2_1 - 1.0) + 
	5.0d-1*(x6_2_1 - 1.0)*(x6_2_1 - 1.0) + 5.0d-1*(y6_2_1 - 1.0)*(y6_2_1 - 1.0) + 
	5.0d-1*(z6_2_1 - 1.0)*(z6_2_1 - 1.0) + 5.0d-1*(x7_2_1 - 1.0)*(x7_2_1 - 1.0) + 
	5.0d-1*(y7_2_1 - 1.0)*(y7_2_1 - 1.0) + 5.0d-1*(z7_2_1 - 1.0)*(z7_2_1 - 1.0) + 
	5.0d-1*(x8_2_1 - 1.0)*(x8_2_1 - 1.0) + 5.0d-1*(y8_2_1 - 1.0)*(y8_2_1 - 1.0) + 
	5.0d-1*(z8_2_1 - 1.0)*(z8_2_1 - 1.0) + 5.0d-1*(x9_2_1 - 1.0)*(x9_2_1 - 1.0) + 
	5.0d-1*(y9_2_1 - 1.0)*(y9_2_1 - 1.0) + 5.0d-1*(z9_2_1 - 1.0)*(z9_2_1 - 1.0) + 
	5.0d-1*(x1_3_1 - 1.0)*(x1_3_1 - 1.0) + 5.0d-1*(y1_3_1 - 1.0)*(y1_3_1 - 1.0) + 
	5.0d-1*(z1_3_1 - 1.0)*(z1_3_1 - 1.0) + 5.0d-1*(x2_3_1 - 1.0)*(x2_3_1 - 1.0) + 
	5.0d-1*(y2_3_1 - 1.0)*(y2_3_1 - 1.0) + 5.0d-1*(z2_3_1 - 1.0)*(z2_3_1 - 1.0) + 
	5.0d-1*(x3_3_1 - 1.0)*(x3_3_1 - 1.0) + 5.0d-1*(y3_3_1 - 1.0)*(y3_3_1 - 1.0) + 
	5.0d-1*(z3_3_1 - 1.0)*(z3_3_1 - 1.0) + 5.0d-1*(x4_3_1 - 1.0)*(x4_3_1 - 1.0) + 
	5.0d-1*(y4_3_1 - 1.0)*(y4_3_1 - 1.0) + 5.0d-1*(z4_3_1 - 1.0)*(z4_3_1 - 1.0) + 
	5.0d-1*(x5_3_1 - 1.0)*(x5_3_1 - 1.0) + 5.0d-1*(y5_3_1 - 1.0)*(y5_3_1 - 1.0) + 
	5.0d-1*(z5_3_1 - 1.0)*(z5_3_1 - 1.0) + 5.0d-1*(x6_3_1 - 1.0)*(x6_3_1 - 1.0) + 
	5.0d-1*(y6_3_1 - 1.0)*(y6_3_1 - 1.0) + 5.0d-1*(z6_3_1 - 1.0)*(z6_3_1 - 1.0) + 
	5.0d-1*(x7_3_1 - 1.0)*(x7_3_1 - 1.0) + 5.0d-1*(y7_3_1 - 1.0)*(y7_3_1 - 1.0) + 
	5.0d-1*(z7_3_1 - 1.0)*(z7_3_1 - 1.0) + 5.0d-1*(x8_3_1 - 1.0)*(x8_3_1 - 1.0) + 
	5.0d-1*(y8_3_1 - 1.0)*(y8_3_1 - 1.0) + 5.0d-1*(z8_3_1 - 1.0)*(z8_3_1 - 1.0) + 
	5.0d-1*(x9_3_1 - 1.0)*(x9_3_1 - 1.0) + 5.0d-1*(y9_3_1 - 1.0)*(y9_3_1 - 1.0) + 
	5.0d-1*(z9_3_1 - 1.0)*(z9_3_1 - 1.0) + 5.0d-1*(x1_4_1 - 1.0)*(x1_4_1 - 1.0) + 
	5.0d-1*(y1_4_1 - 1.0)*(y1_4_1 - 1.0) + 5.0d-1*(z1_4_1 - 1.0)*(z1_4_1 - 1.0) + 
	5.0d-1*(x2_4_1 - 1.0)*(x2_4_1 - 1.0) + 5.0d-1*(y2_4_1 - 1.0)*(y2_4_1 - 1.0) + 
	5.0d-1*(z2_4_1 - 1.0)*(z2_4_1 - 1.0) + 5.0d-1*(x3_4_1 - 1.0)*(x3_4_1 - 1.0) + 
	5.0d-1*(y3_4_1 - 1.0)*(y3_4_1 - 1.0) + 5.0d-1*(z3_4_1 - 1.0)*(z3_4_1 - 1.0) + 
	5.0d-1*(x4_4_1 - 1.0)*(x4_4_1 - 1.0) + 5.0d-1*(y4_4_1 - 1.0)*(y4_4_1 - 1.0) + 
	5.0d-1*(z4_4_1 - 1.0)*(z4_4_1 - 1.0) + 5.0d-1*(x5_4_1 - 1.0)*(x5_4_1 - 1.0) + 
	5.0d-1*(y5_4_1 - 1.0)*(y5_4_1 - 1.0) + 5.0d-1*(z5_4_1 - 1.0)*(z5_4_1 - 1.0) + 
	5.0d-1*(x6_4_1 - 1.0)*(x6_4_1 - 1.0) + 5.0d-1*(y6_4_1 - 1.0)*(y6_4_1 - 1.0) + 
	5.0d-1*(z6_4_1 - 1.0)*(z6_4_1 - 1.0) + 5.0d-1*(x7_4_1 - 1.0)*(x7_4_1 - 1.0) + 
	5.0d-1*(y7_4_1 - 1.0)*(y7_4_1 - 1.0) + 5.0d-1*(z7_4_1 - 1.0)*(z7_4_1 - 1.0) + 
	5.0d-1*(x8_4_1 - 1.0)*(x8_4_1 - 1.0) + 5.0d-1*(y8_4_1 - 1.0)*(y8_4_1 - 1.0) + 
	5.0d-1*(z8_4_1 - 1.0)*(z8_4_1 - 1.0) + 5.0d-1*(x9_4_1 - 1.0)*(x9_4_1 - 1.0) + 
	5.0d-1*(y9_4_1 - 1.0)*(y9_4_1 - 1.0) + 5.0d-1*(z9_4_1 - 1.0)*(z9_4_1 - 1.0) + 
	5.0d-1*(x1_5_1 - 1.0)*(x1_5_1 - 1.0) + 5.0d-1*(y1_5_1 - 1.0)*(y1_5_1 - 1.0) + 
	5.0d-1*(z1_5_1 - 1.0)*(z1_5_1 - 1.0) + 5.0d-1*(x2_5_1 - 1.0)*(x2_5_1 - 1.0) + 
	5.0d-1*(y2_5_1 - 1.0)*(y2_5_1 - 1.0) + 5.0d-1*(z2_5_1 - 1.0)*(z2_5_1 - 1.0) + 
	5.0d-1*(x3_5_1 - 1.0)*(x3_5_1 - 1.0) + 5.0d-1*(y3_5_1 - 1.0)*(y3_5_1 - 1.0) + 
	5.0d-1*(z3_5_1 - 1.0)*(z3_5_1 - 1.0) + 5.0d-1*(x4_5_1 - 1.0)*(x4_5_1 - 1.0) + 
	5.0d-1*(y4_5_1 - 1.0)*(y4_5_1 - 1.0) + 5.0d-1*(z4_5_1 - 1.0)*(z4_5_1 - 1.0) + 
	5.0d-1*(x5_5_1 - 1.0)*(x5_5_1 - 1.0) + 5.0d-1*(y5_5_1 - 1.0)*(y5_5_1 - 1.0) + 
	5.0d-1*(z5_5_1 - 1.0)*(z5_5_1 - 1.0) + 5.0d-1*(x6_5_1 - 1.0)*(x6_5_1 - 1.0) + 
	5.0d-1*(y6_5_1 - 1.0)*(y6_5_1 - 1.0) + 5.0d-1*(z6_5_1 - 1.0)*(z6_5_1 - 1.0) + 
	5.0d-1*(x7_5_1 - 1.0)*(x7_5_1 - 1.0) + 5.0d-1*(y7_5_1 - 1.0)*(y7_5_1 - 1.0) + 
	5.0d-1*(z7_5_1 - 1.0)*(z7_5_1 - 1.0) + 5.0d-1*(x8_5_1 - 1.0)*(x8_5_1 - 1.0) + 
	5.0d-1*(y8_5_1 - 1.0)*(y8_5_1 - 1.0) + 5.0d-1*(z8_5_1 - 1.0)*(z8_5_1 - 1.0) + 
	5.0d-1*(x9_5_1 - 1.0)*(x9_5_1 - 1.0) + 5.0d-1*(y9_5_1 - 1.0)*(y9_5_1 - 1.0) + 
	5.0d-1*(z9_5_1 - 1.0)*(z9_5_1 - 1.0) + 5.0d-1*(x1_6_1 - 1.0)*(x1_6_1 - 1.0) + 
	5.0d-1*(y1_6_1 - 1.0)*(y1_6_1 - 1.0) + 5.0d-1*(z1_6_1 - 1.0)*(z1_6_1 - 1.0) + 
	5.0d-1*(x2_6_1 - 1.0)*(x2_6_1 - 1.0) + 5.0d-1*(y2_6_1 - 1.0)*(y2_6_1 - 1.0) + 
	5.0d-1*(z2_6_1 - 1.0)*(z2_6_1 - 1.0) + 5.0d-1*(x3_6_1 - 1.0)*(x3_6_1 - 1.0) + 
	5.0d-1*(y3_6_1 - 1.0)*(y3_6_1 - 1.0) + 5.0d-1*(z3_6_1 - 1.0)*(z3_6_1 - 1.0) + 
	5.0d-1*(x4_6_1 - 1.0)*(x4_6_1 - 1.0) + 5.0d-1*(y4_6_1 - 1.0)*(y4_6_1 - 1.0) + 
	5.0d-1*(z4_6_1 - 1.0)*(z4_6_1 - 1.0) + 5.0d-1*(x5_6_1 - 1.0)*(x5_6_1 - 1.0) + 
	5.0d-1*(y5_6_1 - 1.0)*(y5_6_1 - 1.0) + 5.0d-1*(z5_6_1 - 1.0)*(z5_6_1 - 1.0) + 
	5.0d-1*(x6_6_1 - 1.0)*(x6_6_1 - 1.0) + 5.0d-1*(y6_6_1 - 1.0)*(y6_6_1 - 1.0) + 
	5.0d-1*(z6_6_1 - 1.0)*(z6_6_1 - 1.0) + 5.0d-1*(x7_6_1 - 1.0)*(x7_6_1 - 1.0) + 
	5.0d-1*(y7_6_1 - 1.0)*(y7_6_1 - 1.0) + 5.0d-1*(z7_6_1 - 1.0)*(z7_6_1 - 1.0) + 
	5.0d-1*(x8_6_1 - 1.0)*(x8_6_1 - 1.0) + 5.0d-1*(y8_6_1 - 1.0)*(y8_6_1 - 1.0) + 
	5.0d-1*(z8_6_1 - 1.0)*(z8_6_1 - 1.0) + 5.0d-1*(x9_6_1 - 1.0)*(x9_6_1 - 1.0) + 
	5.0d-1*(y9_6_1 - 1.0)*(y9_6_1 - 1.0) + 5.0d-1*(z9_6_1 - 1.0)*(z9_6_1 - 1.0) + 
	5.0d-1*(x1_7_1 - 1.0)*(x1_7_1 - 1.0) + 5.0d-1*(y1_7_1 - 1.0)*(y1_7_1 - 1.0) + 
	5.0d-1*(z1_7_1 - 1.0)*(z1_7_1 - 1.0) + 5.0d-1*(x2_7_1 - 1.0)*(x2_7_1 - 1.0) + 
	5.0d-1*(y2_7_1 - 1.0)*(y2_7_1 - 1.0) + 5.0d-1*(z2_7_1 - 1.0)*(z2_7_1 - 1.0) + 
	5.0d-1*(x3_7_1 - 1.0)*(x3_7_1 - 1.0) + 5.0d-1*(y3_7_1 - 1.0)*(y3_7_1 - 1.0) + 
	5.0d-1*(z3_7_1 - 1.0)*(z3_7_1 - 1.0) + 5.0d-1*(x4_7_1 - 1.0)*(x4_7_1 - 1.0) + 
	5.0d-1*(y4_7_1 - 1.0)*(y4_7_1 - 1.0) + 5.0d-1*(z4_7_1 - 1.0)*(z4_7_1 - 1.0) + 
	5.0d-1*(x5_7_1 - 1.0)*(x5_7_1 - 1.0) + 5.0d-1*(y5_7_1 - 1.0)*(y5_7_1 - 1.0) + 
	5.0d-1*(z5_7_1 - 1.0)*(z5_7_1 - 1.0) + 5.0d-1*(x6_7_1 - 1.0)*(x6_7_1 - 1.0) + 
	5.0d-1*(y6_7_1 - 1.0)*(y6_7_1 - 1.0) + 5.0d-1*(z6_7_1 - 1.0)*(z6_7_1 - 1.0) + 
	5.0d-1*(x7_7_1 - 1.0)*(x7_7_1 - 1.0) + 5.0d-1*(y7_7_1 - 1.0)*(y7_7_1 - 1.0) + 
	5.0d-1*(z7_7_1 - 1.0)*(z7_7_1 - 1.0) + 5.0d-1*(x8_7_1 - 1.0)*(x8_7_1 - 1.0) + 
	5.0d-1*(y8_7_1 - 1.0)*(y8_7_1 - 1.0) + 5.0d-1*(z8_7_1 - 1.0)*(z8_7_1 - 1.0) + 
	5.0d-1*(x9_7_1 - 1.0)*(x9_7_1 - 1.0) + 5.0d-1*(y9_7_1 - 1.0)*(y9_7_1 - 1.0) + 
	5.0d-1*(z9_7_1 - 1.0)*(z9_7_1 - 1.0) + 5.0d-1*(x1_8_1 - 1.0)*(x1_8_1 - 1.0) + 
	5.0d-1*(y1_8_1 - 1.0)*(y1_8_1 - 1.0) + 5.0d-1*(z1_8_1 - 1.0)*(z1_8_1 - 1.0) + 
	5.0d-1*(x2_8_1 - 1.0)*(x2_8_1 - 1.0) + 5.0d-1*(y2_8_1 - 1.0)*(y2_8_1 - 1.0) + 
	5.0d-1*(z2_8_1 - 1.0)*(z2_8_1 - 1.0) + 5.0d-1*(x3_8_1 - 1.0)*(x3_8_1 - 1.0) + 
	5.0d-1*(y3_8_1 - 1.0)*(y3_8_1 - 1.0) + 5.0d-1*(z3_8_1 - 1.0)*(z3_8_1 - 1.0) + 
	5.0d-1*(x4_8_1 - 1.0)*(x4_8_1 - 1.0) + 5.0d-1*(y4_8_1 - 1.0)*(y4_8_1 - 1.0) + 
	5.0d-1*(z4_8_1 - 1.0)*(z4_8_1 - 1.0) + 5.0d-1*(x5_8_1 - 1.0)*(x5_8_1 - 1.0) + 
	5.0d-1*(y5_8_1 - 1.0)*(y5_8_1 - 1.0) + 5.0d-1*(z5_8_1 - 1.0)*(z5_8_1 - 1.0) + 
	5.0d-1*(x6_8_1 - 1.0)*(x6_8_1 - 1.0) + 5.0d-1*(y6_8_1 - 1.0)*(y6_8_1 - 1.0) + 
	5.0d-1*(z6_8_1 - 1.0)*(z6_8_1 - 1.0) + 5.0d-1*(x7_8_1 - 1.0)*(x7_8_1 - 1.0) + 
	5.0d-1*(y7_8_1 - 1.0)*(y7_8_1 - 1.0) + 5.0d-1*(z7_8_1 - 1.0)*(z7_8_1 - 1.0) + 
	5.0d-1*(x8_8_1 - 1.0)*(x8_8_1 - 1.0) + 5.0d-1*(y8_8_1 - 1.0)*(y8_8_1 - 1.0) + 
	5.0d-1*(z8_8_1 - 1.0)*(z8_8_1 - 1.0) + 5.0d-1*(x9_8_1 - 1.0)*(x9_8_1 - 1.0) + 
	5.0d-1*(y9_8_1 - 1.0)*(y9_8_1 - 1.0) + 5.0d-1*(z9_8_1 - 1.0)*(z9_8_1 - 1.0) + 
	5.0d-1*(x1_9_1 - 1.0)*(x1_9_1 - 1.0) + 5.0d-1*(y1_9_1 - 1.0)*(y1_9_1 - 1.0) + 
	5.0d-1*(z1_9_1 - 1.0)*(z1_9_1 - 1.0) + 5.0d-1*(x2_9_1 - 1.0)*(x2_9_1 - 1.0) + 
	5.0d-1*(y2_9_1 - 1.0)*(y2_9_1 - 1.0) + 5.0d-1*(z2_9_1 - 1.0)*(z2_9_1 - 1.0) + 
	5.0d-1*(x3_9_1 - 1.0)*(x3_9_1 - 1.0) + 5.0d-1*(y3_9_1 - 1.0)*(y3_9_1 - 1.0) + 
	5.0d-1*(z3_9_1 - 1.0)*(z3_9_1 - 1.0) + 5.0d-1*(x4_9_1 - 1.0)*(x4_9_1 - 1.0) + 
	5.0d-1*(y4_9_1 - 1.0)*(y4_9_1 - 1.0) + 5.0d-1*(z4_9_1 - 1.0)*(z4_9_1 - 1.0) + 
	5.0d-1*(x5_9_1 - 1.0)*(x5_9_1 - 1.0) + 5.0d-1*(y5_9_1 - 1.0)*(y5_9_1 - 1.0) + 
	5.0d-1*(z5_9_1 - 1.0)*(z5_9_1 - 1.0) + 5.0d-1*(x6_9_1 - 1.0)*(x6_9_1 - 1.0) + 
	5.0d-1*(y6_9_1 - 1.0)*(y6_9_1 - 1.0) + 5.0d-1*(z6_9_1 - 1.0)*(z6_9_1 - 1.0) + 
	5.0d-1*(x7_9_1 - 1.0)*(x7_9_1 - 1.0) + 5.0d-1*(y7_9_1 - 1.0)*(y7_9_1 - 1.0) + 
	5.0d-1*(z7_9_1 - 1.0)*(z7_9_1 - 1.0) + 5.0d-1*(x8_9_1 - 1.0)*(x8_9_1 - 1.0) + 
	5.0d-1*(y8_9_1 - 1.0)*(y8_9_1 - 1.0) + 5.0d-1*(z8_9_1 - 1.0)*(z8_9_1 - 1.0) + 
	5.0d-1*(x9_9_1 - 1.0)*(x9_9_1 - 1.0) + 5.0d-1*(y9_9_1 - 1.0)*(y9_9_1 - 1.0) + 
	5.0d-1*(z9_9_1 - 1.0)*(z9_9_1 - 1.0) + 5.0d-1*(x1_1_2 - 1.0)*(x1_1_2 - 1.0) + 
	5.0d-1*(y1_1_2 - 1.0)*(y1_1_2 - 1.0) + 5.0d-1*(z1_1_2 - 1.0)*(z1_1_2 - 1.0) + 
	5.0d-1*(x2_1_2 - 1.0)*(x2_1_2 - 1.0) + 5.0d-1*(y2_1_2 - 1.0)*(y2_1_2 - 1.0) + 
	5.0d-1*(z2_1_2 - 1.0)*(z2_1_2 - 1.0) + 5.0d-1*(x3_1_2 - 1.0)*(x3_1_2 - 1.0) + 
	5.0d-1*(y3_1_2 - 1.0)*(y3_1_2 - 1.0) + 5.0d-1*(z3_1_2 - 1.0)*(z3_1_2 - 1.0) + 
	5.0d-1*(x4_1_2 - 1.0)*(x4_1_2 - 1.0) + 5.0d-1*(y4_1_2 - 1.0)*(y4_1_2 - 1.0) + 
	5.0d-1*(z4_1_2 - 1.0)*(z4_1_2 - 1.0) + 5.0d-1*(x5_1_2 - 1.0)*(x5_1_2 - 1.0) + 
	5.0d-1*(y5_1_2 - 1.0)*(y5_1_2 - 1.0) + 5.0d-1*(z5_1_2 - 1.0)*(z5_1_2 - 1.0) + 
	5.0d-1*(x6_1_2 - 1.0)*(x6_1_2 - 1.0) + 5.0d-1*(y6_1_2 - 1.0)*(y6_1_2 - 1.0) + 
	5.0d-1*(z6_1_2 - 1.0)*(z6_1_2 - 1.0) + 5.0d-1*(x7_1_2 - 1.0)*(x7_1_2 - 1.0) + 
	5.0d-1*(y7_1_2 - 1.0)*(y7_1_2 - 1.0) + 5.0d-1*(z7_1_2 - 1.0)*(z7_1_2 - 1.0) + 
	5.0d-1*(x8_1_2 - 1.0)*(x8_1_2 - 1.0) + 5.0d-1*(y8_1_2 - 1.0)*(y8_1_2 - 1.0) + 
	5.0d-1*(z8_1_2 - 1.0)*(z8_1_2 - 1.0) + 5.0d-1*(x9_1_2 - 1.0)*(x9_1_2 - 1.0) + 
	5.0d-1*(y9_1_2 - 1.0)*(y9_1_2 - 1.0) + 5.0d-1*(z9_1_2 - 1.0)*(z9_1_2 - 1.0) + 
	5.0d-1*(x1_2_2 - 1.0)*(x1_2_2 - 1.0) + 5.0d-1*(y1_2_2 - 1.0)*(y1_2_2 - 1.0) + 
	5.0d-1*(z1_2_2 - 1.0)*(z1_2_2 - 1.0) + 5.0d-1*(x2_2_2 - 1.0)*(x2_2_2 - 1.0) + 
	5.0d-1*(y2_2_2 - 1.0)*(y2_2_2 - 1.0) + 5.0d-1*(z2_2_2 - 1.0)*(z2_2_2 - 1.0) + 
	5.0d-1*(x3_2_2 - 1.0)*(x3_2_2 - 1.0) + 5.0d-1*(y3_2_2 - 1.0)*(y3_2_2 - 1.0) + 
	5.0d-1*(z3_2_2 - 1.0)*(z3_2_2 - 1.0) + 5.0d-1*(x4_2_2 - 1.0)*(x4_2_2 - 1.0) + 
	5.0d-1*(y4_2_2 - 1.0)*(y4_2_2 - 1.0) + 5.0d-1*(z4_2_2 - 1.0)*(z4_2_2 - 1.0) + 
	5.0d-1*(x5_2_2 - 1.0)*(x5_2_2 - 1.0) + 5.0d-1*(y5_2_2 - 1.0)*(y5_2_2 - 1.0) + 
	5.0d-1*(z5_2_2 - 1.0)*(z5_2_2 - 1.0) + 5.0d-1*(x6_2_2 - 1.0)*(x6_2_2 - 1.0) + 
	5.0d-1*(y6_2_2 - 1.0)*(y6_2_2 - 1.0) + 5.0d-1*(z6_2_2 - 1.0)*(z6_2_2 - 1.0) + 
	5.0d-1*(x7_2_2 - 1.0)*(x7_2_2 - 1.0) + 5.0d-1*(y7_2_2 - 1.0)*(y7_2_2 - 1.0) + 
	5.0d-1*(z7_2_2 - 1.0)*(z7_2_2 - 1.0) + 5.0d-1*(x8_2_2 - 1.0)*(x8_2_2 - 1.0) + 
	5.0d-1*(y8_2_2 - 1.0)*(y8_2_2 - 1.0) + 5.0d-1*(z8_2_2 - 1.0)*(z8_2_2 - 1.0) + 
	5.0d-1*(x9_2_2 - 1.0)*(x9_2_2 - 1.0) + 5.0d-1*(y9_2_2 - 1.0)*(y9_2_2 - 1.0) + 
	5.0d-1*(z9_2_2 - 1.0)*(z9_2_2 - 1.0) + 5.0d-1*(x1_3_2 - 1.0)*(x1_3_2 - 1.0) + 
	5.0d-1*(y1_3_2 - 1.0)*(y1_3_2 - 1.0) + 5.0d-1*(z1_3_2 - 1.0)*(z1_3_2 - 1.0) + 
	5.0d-1*(x2_3_2 - 1.0)*(x2_3_2 - 1.0) + 5.0d-1*(y2_3_2 - 1.0)*(y2_3_2 - 1.0) + 
	5.0d-1*(z2_3_2 - 1.0)*(z2_3_2 - 1.0) + 5.0d-1*(x3_3_2 - 1.0)*(x3_3_2 - 1.0) + 
	5.0d-1*(y3_3_2 - 1.0)*(y3_3_2 - 1.0) + 5.0d-1*(z3_3_2 - 1.0)*(z3_3_2 - 1.0) + 
	5.0d-1*(x4_3_2 - 1.0)*(x4_3_2 - 1.0) + 5.0d-1*(y4_3_2 - 1.0)*(y4_3_2 - 1.0) + 
	5.0d-1*(z4_3_2 - 1.0)*(z4_3_2 - 1.0) + 5.0d-1*(x5_3_2 - 1.0)*(x5_3_2 - 1.0) + 
	5.0d-1*(y5_3_2 - 1.0)*(y5_3_2 - 1.0) + 5.0d-1*(z5_3_2 - 1.0)*(z5_3_2 - 1.0) + 
	5.0d-1*(x6_3_2 - 1.0)*(x6_3_2 - 1.0) + 5.0d-1*(y6_3_2 - 1.0)*(y6_3_2 - 1.0) + 
	5.0d-1*(z6_3_2 - 1.0)*(z6_3_2 - 1.0) + 5.0d-1*(x7_3_2 - 1.0)*(x7_3_2 - 1.0) + 
	5.0d-1*(y7_3_2 - 1.0)*(y7_3_2 - 1.0) + 5.0d-1*(z7_3_2 - 1.0)*(z7_3_2 - 1.0) + 
	5.0d-1*(x8_3_2 - 1.0)*(x8_3_2 - 1.0) + 5.0d-1*(y8_3_2 - 1.0)*(y8_3_2 - 1.0) + 
	5.0d-1*(z8_3_2 - 1.0)*(z8_3_2 - 1.0) + 5.0d-1*(x9_3_2 - 1.0)*(x9_3_2 - 1.0) + 
	5.0d-1*(y9_3_2 - 1.0)*(y9_3_2 - 1.0) + 5.0d-1*(z9_3_2 - 1.0)*(z9_3_2 - 1.0) + 
	5.0d-1*(x1_4_2 - 1.0)*(x1_4_2 - 1.0) + 5.0d-1*(y1_4_2 - 1.0)*(y1_4_2 - 1.0) + 
	5.0d-1*(z1_4_2 - 1.0)*(z1_4_2 - 1.0) + 5.0d-1*(x2_4_2 - 1.0)*(x2_4_2 - 1.0) + 
	5.0d-1*(y2_4_2 - 1.0)*(y2_4_2 - 1.0) + 5.0d-1*(z2_4_2 - 1.0)*(z2_4_2 - 1.0) + 
	5.0d-1*(x3_4_2 - 1.0)*(x3_4_2 - 1.0) + 5.0d-1*(y3_4_2 - 1.0)*(y3_4_2 - 1.0) + 
	5.0d-1*(z3_4_2 - 1.0)*(z3_4_2 - 1.0) + 5.0d-1*(x4_4_2 - 1.0)*(x4_4_2 - 1.0) + 
	5.0d-1*(y4_4_2 - 1.0)*(y4_4_2 - 1.0) + 5.0d-1*(z4_4_2 - 1.0)*(z4_4_2 - 1.0) + 
	5.0d-1*(x5_4_2 - 1.0)*(x5_4_2 - 1.0) + 5.0d-1*(y5_4_2 - 1.0)*(y5_4_2 - 1.0) + 
	5.0d-1*(z5_4_2 - 1.0)*(z5_4_2 - 1.0) + 5.0d-1*(x6_4_2 - 1.0)*(x6_4_2 - 1.0) + 
	5.0d-1*(y6_4_2 - 1.0)*(y6_4_2 - 1.0) + 5.0d-1*(z6_4_2 - 1.0)*(z6_4_2 - 1.0) + 
	5.0d-1*(x7_4_2 - 1.0)*(x7_4_2 - 1.0) + 5.0d-1*(y7_4_2 - 1.0)*(y7_4_2 - 1.0) + 
	5.0d-1*(z7_4_2 - 1.0)*(z7_4_2 - 1.0) + 5.0d-1*(x8_4_2 - 1.0)*(x8_4_2 - 1.0) + 
	5.0d-1*(y8_4_2 - 1.0)*(y8_4_2 - 1.0) + 5.0d-1*(z8_4_2 - 1.0)*(z8_4_2 - 1.0) + 
	5.0d-1*(x9_4_2 - 1.0)*(x9_4_2 - 1.0) + 5.0d-1*(y9_4_2 - 1.0)*(y9_4_2 - 1.0) + 
	5.0d-1*(z9_4_2 - 1.0)*(z9_4_2 - 1.0) + 5.0d-1*(x1_5_2 - 1.0)*(x1_5_2 - 1.0) + 
	5.0d-1*(y1_5_2 - 1.0)*(y1_5_2 - 1.0) + 5.0d-1*(z1_5_2 - 1.0)*(z1_5_2 - 1.0) + 
	5.0d-1*(x2_5_2 - 1.0)*(x2_5_2 - 1.0) + 5.0d-1*(y2_5_2 - 1.0)*(y2_5_2 - 1.0) + 
	5.0d-1*(z2_5_2 - 1.0)*(z2_5_2 - 1.0) + 5.0d-1*(x3_5_2 - 1.0)*(x3_5_2 - 1.0) + 
	5.0d-1*(y3_5_2 - 1.0)*(y3_5_2 - 1.0) + 5.0d-1*(z3_5_2 - 1.0)*(z3_5_2 - 1.0) + 
	5.0d-1*(x4_5_2 - 1.0)*(x4_5_2 - 1.0) + 5.0d-1*(y4_5_2 - 1.0)*(y4_5_2 - 1.0) + 
	5.0d-1*(z4_5_2 - 1.0)*(z4_5_2 - 1.0) + 5.0d-1*(x5_5_2 - 1.0)*(x5_5_2 - 1.0) + 
	5.0d-1*(y5_5_2 - 1.0)*(y5_5_2 - 1.0) + 5.0d-1*(z5_5_2 - 1.0)*(z5_5_2 - 1.0) + 
	5.0d-1*(x6_5_2 - 1.0)*(x6_5_2 - 1.0) + 5.0d-1*(y6_5_2 - 1.0)*(y6_5_2 - 1.0) + 
	5.0d-1*(z6_5_2 - 1.0)*(z6_5_2 - 1.0) + 5.0d-1*(x7_5_2 - 1.0)*(x7_5_2 - 1.0) + 
	5.0d-1*(y7_5_2 - 1.0)*(y7_5_2 - 1.0) + 5.0d-1*(z7_5_2 - 1.0)*(z7_5_2 - 1.0) + 
	5.0d-1*(x8_5_2 - 1.0)*(x8_5_2 - 1.0) + 5.0d-1*(y8_5_2 - 1.0)*(y8_5_2 - 1.0) + 
	5.0d-1*(z8_5_2 - 1.0)*(z8_5_2 - 1.0) + 5.0d-1*(x9_5_2 - 1.0)*(x9_5_2 - 1.0) + 
	5.0d-1*(y9_5_2 - 1.0)*(y9_5_2 - 1.0) + 5.0d-1*(z9_5_2 - 1.0)*(z9_5_2 - 1.0) + 
	5.0d-1*(x1_6_2 - 1.0)*(x1_6_2 - 1.0) + 5.0d-1*(y1_6_2 - 1.0)*(y1_6_2 - 1.0) + 
	5.0d-1*(z1_6_2 - 1.0)*(z1_6_2 - 1.0) + 5.0d-1*(x2_6_2 - 1.0)*(x2_6_2 - 1.0) + 
	5.0d-1*(y2_6_2 - 1.0)*(y2_6_2 - 1.0) + 5.0d-1*(z2_6_2 - 1.0)*(z2_6_2 - 1.0) + 
	5.0d-1*(x3_6_2 - 1.0)*(x3_6_2 - 1.0) + 5.0d-1*(y3_6_2 - 1.0)*(y3_6_2 - 1.0) + 
	5.0d-1*(z3_6_2 - 1.0)*(z3_6_2 - 1.0) + 5.0d-1*(x4_6_2 - 1.0)*(x4_6_2 - 1.0) + 
	5.0d-1*(y4_6_2 - 1.0)*(y4_6_2 - 1.0) + 5.0d-1*(z4_6_2 - 1.0)*(z4_6_2 - 1.0) + 
	5.0d-1*(x5_6_2 - 1.0)*(x5_6_2 - 1.0) + 5.0d-1*(y5_6_2 - 1.0)*(y5_6_2 - 1.0) + 
	5.0d-1*(z5_6_2 - 1.0)*(z5_6_2 - 1.0) + 5.0d-1*(x6_6_2 - 1.0)*(x6_6_2 - 1.0) + 
	5.0d-1*(y6_6_2 - 1.0)*(y6_6_2 - 1.0) + 5.0d-1*(z6_6_2 - 1.0)*(z6_6_2 - 1.0) + 
	5.0d-1*(x7_6_2 - 1.0)*(x7_6_2 - 1.0) + 5.0d-1*(y7_6_2 - 1.0)*(y7_6_2 - 1.0) + 
	5.0d-1*(z7_6_2 - 1.0)*(z7_6_2 - 1.0) + 5.0d-1*(x8_6_2 - 1.0)*(x8_6_2 - 1.0) + 
	5.0d-1*(y8_6_2 - 1.0)*(y8_6_2 - 1.0) + 5.0d-1*(z8_6_2 - 1.0)*(z8_6_2 - 1.0) + 
	5.0d-1*(x9_6_2 - 1.0)*(x9_6_2 - 1.0) + 5.0d-1*(y9_6_2 - 1.0)*(y9_6_2 - 1.0) + 
	5.0d-1*(z9_6_2 - 1.0)*(z9_6_2 - 1.0) + 5.0d-1*(x1_7_2 - 1.0)*(x1_7_2 - 1.0) + 
	5.0d-1*(y1_7_2 - 1.0)*(y1_7_2 - 1.0) + 5.0d-1*(z1_7_2 - 1.0)*(z1_7_2 - 1.0) + 
	5.0d-1*(x2_7_2 - 1.0)*(x2_7_2 - 1.0) + 5.0d-1*(y2_7_2 - 1.0)*(y2_7_2 - 1.0) + 
	5.0d-1*(z2_7_2 - 1.0)*(z2_7_2 - 1.0) + 5.0d-1*(x3_7_2 - 1.0)*(x3_7_2 - 1.0) + 
	5.0d-1*(y3_7_2 - 1.0)*(y3_7_2 - 1.0) + 5.0d-1*(z3_7_2 - 1.0)*(z3_7_2 - 1.0) + 
	5.0d-1*(x4_7_2 - 1.0)*(x4_7_2 - 1.0) + 5.0d-1*(y4_7_2 - 1.0)*(y4_7_2 - 1.0) + 
	5.0d-1*(z4_7_2 - 1.0)*(z4_7_2 - 1.0) + 5.0d-1*(x5_7_2 - 1.0)*(x5_7_2 - 1.0) + 
	5.0d-1*(y5_7_2 - 1.0)*(y5_7_2 - 1.0) + 5.0d-1*(z5_7_2 - 1.0)*(z5_7_2 - 1.0) + 
	5.0d-1*(x6_7_2 - 1.0)*(x6_7_2 - 1.0) + 5.0d-1*(y6_7_2 - 1.0)*(y6_7_2 - 1.0) + 
	5.0d-1*(z6_7_2 - 1.0)*(z6_7_2 - 1.0) + 5.0d-1*(x7_7_2 - 1.0)*(x7_7_2 - 1.0) + 
	5.0d-1*(y7_7_2 - 1.0)*(y7_7_2 - 1.0) + 5.0d-1*(z7_7_2 - 1.0)*(z7_7_2 - 1.0) + 
	5.0d-1*(x8_7_2 - 1.0)*(x8_7_2 - 1.0) + 5.0d-1*(y8_7_2 - 1.0)*(y8_7_2 - 1.0) + 
	5.0d-1*(z8_7_2 - 1.0)*(z8_7_2 - 1.0) + 5.0d-1*(x9_7_2 - 1.0)*(x9_7_2 - 1.0) + 
	5.0d-1*(y9_7_2 - 1.0)*(y9_7_2 - 1.0) + 5.0d-1*(z9_7_2 - 1.0)*(z9_7_2 - 1.0) + 
	5.0d-1*(x1_8_2 - 1.0)*(x1_8_2 - 1.0) + 5.0d-1*(y1_8_2 - 1.0)*(y1_8_2 - 1.0) + 
	5.0d-1*(z1_8_2 - 1.0)*(z1_8_2 - 1.0) + 5.0d-1*(x2_8_2 - 1.0)*(x2_8_2 - 1.0) + 
	5.0d-1*(y2_8_2 - 1.0)*(y2_8_2 - 1.0) + 5.0d-1*(z2_8_2 - 1.0)*(z2_8_2 - 1.0) + 
	5.0d-1*(x3_8_2 - 1.0)*(x3_8_2 - 1.0) + 5.0d-1*(y3_8_2 - 1.0)*(y3_8_2 - 1.0) + 
	5.0d-1*(z3_8_2 - 1.0)*(z3_8_2 - 1.0) + 5.0d-1*(x4_8_2 - 1.0)*(x4_8_2 - 1.0) + 
	5.0d-1*(y4_8_2 - 1.0)*(y4_8_2 - 1.0) + 5.0d-1*(z4_8_2 - 1.0)*(z4_8_2 - 1.0) + 
	5.0d-1*(x5_8_2 - 1.0)*(x5_8_2 - 1.0) + 5.0d-1*(y5_8_2 - 1.0)*(y5_8_2 - 1.0) + 
	5.0d-1*(z5_8_2 - 1.0)*(z5_8_2 - 1.0) + 5.0d-1*(x6_8_2 - 1.0)*(x6_8_2 - 1.0) + 
	5.0d-1*(y6_8_2 - 1.0)*(y6_8_2 - 1.0) + 5.0d-1*(z6_8_2 - 1.0)*(z6_8_2 - 1.0) + 
	5.0d-1*(x7_8_2 - 1.0)*(x7_8_2 - 1.0) + 5.0d-1*(y7_8_2 - 1.0)*(y7_8_2 - 1.0) + 
	5.0d-1*(z7_8_2 - 1.0)*(z7_8_2 - 1.0) + 5.0d-1*(x8_8_2 - 1.0)*(x8_8_2 - 1.0) + 
	5.0d-1*(y8_8_2 - 1.0)*(y8_8_2 - 1.0) + 5.0d-1*(z8_8_2 - 1.0)*(z8_8_2 - 1.0) + 
	5.0d-1*(x9_8_2 - 1.0)*(x9_8_2 - 1.0) + 5.0d-1*(y9_8_2 - 1.0)*(y9_8_2 - 1.0) + 
	5.0d-1*(z9_8_2 - 1.0)*(z9_8_2 - 1.0) + 5.0d-1*(x1_9_2 - 1.0)*(x1_9_2 - 1.0) + 
	5.0d-1*(y1_9_2 - 1.0)*(y1_9_2 - 1.0) + 5.0d-1*(z1_9_2 - 1.0)*(z1_9_2 - 1.0) + 
	5.0d-1*(x2_9_2 - 1.0)*(x2_9_2 - 1.0) + 5.0d-1*(y2_9_2 - 1.0)*(y2_9_2 - 1.0) + 
	5.0d-1*(z2_9_2 - 1.0)*(z2_9_2 - 1.0) + 5.0d-1*(x3_9_2 - 1.0)*(x3_9_2 - 1.0) + 
	5.0d-1*(y3_9_2 - 1.0)*(y3_9_2 - 1.0) + 5.0d-1*(z3_9_2 - 1.0)*(z3_9_2 - 1.0) + 
	5.0d-1*(x4_9_2 - 1.0)*(x4_9_2 - 1.0) + 5.0d-1*(y4_9_2 - 1.0)*(y4_9_2 - 1.0) + 
	5.0d-1*(z4_9_2 - 1.0)*(z4_9_2 - 1.0) + 5.0d-1*(x5_9_2 - 1.0)*(x5_9_2 - 1.0) + 
	5.0d-1*(y5_9_2 - 1.0)*(y5_9_2 - 1.0) + 5.0d-1*(z5_9_2 - 1.0)*(z5_9_2 - 1.0) + 
	5.0d-1*(x6_9_2 - 1.0)*(x6_9_2 - 1.0) + 5.0d-1*(y6_9_2 - 1.0)*(y6_9_2 - 1.0) + 
	5.0d-1*(z6_9_2 - 1.0)*(z6_9_2 - 1.0) + 5.0d-1*(x7_9_2 - 1.0)*(x7_9_2 - 1.0) + 
	5.0d-1*(y7_9_2 - 1.0)*(y7_9_2 - 1.0) + 5.0d-1*(z7_9_2 - 1.0)*(z7_9_2 - 1.0) + 
	5.0d-1*(x8_9_2 - 1.0)*(x8_9_2 - 1.0) + 5.0d-1*(y8_9_2 - 1.0)*(y8_9_2 - 1.0) + 
	5.0d-1*(z8_9_2 - 1.0)*(z8_9_2 - 1.0) + 5.0d-1*(x9_9_2 - 1.0)*(x9_9_2 - 1.0) + 
	5.0d-1*(y9_9_2 - 1.0)*(y9_9_2 - 1.0) + 5.0d-1*(z9_9_2 - 1.0)*(z9_9_2 - 1.0) + 
	5.0d-1*(x1_1_3 - 1.0)*(x1_1_3 - 1.0) + 5.0d-1*(y1_1_3 - 1.0)*(y1_1_3 - 1.0) + 
	5.0d-1*(z1_1_3 - 1.0)*(z1_1_3 - 1.0) + 5.0d-1*(x2_1_3 - 1.0)*(x2_1_3 - 1.0) + 
	5.0d-1*(y2_1_3 - 1.0)*(y2_1_3 - 1.0) + 5.0d-1*(z2_1_3 - 1.0)*(z2_1_3 - 1.0) + 
	5.0d-1*(x3_1_3 - 1.0)*(x3_1_3 - 1.0) + 5.0d-1*(y3_1_3 - 1.0)*(y3_1_3 - 1.0) + 
	5.0d-1*(z3_1_3 - 1.0)*(z3_1_3 - 1.0) + 5.0d-1*(x4_1_3 - 1.0)*(x4_1_3 - 1.0) + 
	5.0d-1*(y4_1_3 - 1.0)*(y4_1_3 - 1.0) + 5.0d-1*(z4_1_3 - 1.0)*(z4_1_3 - 1.0) + 
	5.0d-1*(x5_1_3 - 1.0)*(x5_1_3 - 1.0) + 5.0d-1*(y5_1_3 - 1.0)*(y5_1_3 - 1.0) + 
	5.0d-1*(z5_1_3 - 1.0)*(z5_1_3 - 1.0) + 5.0d-1*(x6_1_3 - 1.0)*(x6_1_3 - 1.0) + 
	5.0d-1*(y6_1_3 - 1.0)*(y6_1_3 - 1.0) + 5.0d-1*(z6_1_3 - 1.0)*(z6_1_3 - 1.0) + 
	5.0d-1*(x7_1_3 - 1.0)*(x7_1_3 - 1.0) + 5.0d-1*(y7_1_3 - 1.0)*(y7_1_3 - 1.0) + 
	5.0d-1*(z7_1_3 - 1.0)*(z7_1_3 - 1.0) + 5.0d-1*(x8_1_3 - 1.0)*(x8_1_3 - 1.0) + 
	5.0d-1*(y8_1_3 - 1.0)*(y8_1_3 - 1.0) + 5.0d-1*(z8_1_3 - 1.0)*(z8_1_3 - 1.0) + 
	5.0d-1*(x9_1_3 - 1.0)*(x9_1_3 - 1.0) + 5.0d-1*(y9_1_3 - 1.0)*(y9_1_3 - 1.0) + 
	5.0d-1*(z9_1_3 - 1.0)*(z9_1_3 - 1.0) + 5.0d-1*(x1_2_3 - 1.0)*(x1_2_3 - 1.0) + 
	5.0d-1*(y1_2_3 - 1.0)*(y1_2_3 - 1.0) + 5.0d-1*(z1_2_3 - 1.0)*(z1_2_3 - 1.0) + 
	5.0d-1*(x2_2_3 - 1.0)*(x2_2_3 - 1.0) + 5.0d-1*(y2_2_3 - 1.0)*(y2_2_3 - 1.0) + 
	5.0d-1*(z2_2_3 - 1.0)*(z2_2_3 - 1.0) + 5.0d-1*(x3_2_3 - 1.0)*(x3_2_3 - 1.0) + 
	5.0d-1*(y3_2_3 - 1.0)*(y3_2_3 - 1.0) + 5.0d-1*(z3_2_3 - 1.0)*(z3_2_3 - 1.0) + 
	5.0d-1*(x4_2_3 - 1.0)*(x4_2_3 - 1.0) + 5.0d-1*(y4_2_3 - 1.0)*(y4_2_3 - 1.0) + 
	5.0d-1*(z4_2_3 - 1.0)*(z4_2_3 - 1.0) + 5.0d-1*(x5_2_3 - 1.0)*(x5_2_3 - 1.0) + 
	5.0d-1*(y5_2_3 - 1.0)*(y5_2_3 - 1.0) + 5.0d-1*(z5_2_3 - 1.0)*(z5_2_3 - 1.0) + 
	5.0d-1*(x6_2_3 - 1.0)*(x6_2_3 - 1.0) + 5.0d-1*(y6_2_3 - 1.0)*(y6_2_3 - 1.0) + 
	5.0d-1*(z6_2_3 - 1.0)*(z6_2_3 - 1.0) + 5.0d-1*(x7_2_3 - 1.0)*(x7_2_3 - 1.0) + 
	5.0d-1*(y7_2_3 - 1.0)*(y7_2_3 - 1.0) + 5.0d-1*(z7_2_3 - 1.0)*(z7_2_3 - 1.0) + 
	5.0d-1*(x8_2_3 - 1.0)*(x8_2_3 - 1.0) + 5.0d-1*(y8_2_3 - 1.0)*(y8_2_3 - 1.0) + 
	5.0d-1*(z8_2_3 - 1.0)*(z8_2_3 - 1.0) + 5.0d-1*(x9_2_3 - 1.0)*(x9_2_3 - 1.0) + 
	5.0d-1*(y9_2_3 - 1.0)*(y9_2_3 - 1.0) + 5.0d-1*(z9_2_3 - 1.0)*(z9_2_3 - 1.0) + 
	5.0d-1*(x1_3_3 - 1.0)*(x1_3_3 - 1.0) + 5.0d-1*(y1_3_3 - 1.0)*(y1_3_3 - 1.0) + 
	5.0d-1*(z1_3_3 - 1.0)*(z1_3_3 - 1.0) + 5.0d-1*(x2_3_3 - 1.0)*(x2_3_3 - 1.0) + 
	5.0d-1*(y2_3_3 - 1.0)*(y2_3_3 - 1.0) + 5.0d-1*(z2_3_3 - 1.0)*(z2_3_3 - 1.0) + 
	5.0d-1*(x3_3_3 - 1.0)*(x3_3_3 - 1.0) + 5.0d-1*(y3_3_3 - 1.0)*(y3_3_3 - 1.0) + 
	5.0d-1*(z3_3_3 - 1.0)*(z3_3_3 - 1.0) + 5.0d-1*(x4_3_3 - 1.0)*(x4_3_3 - 1.0) + 
	5.0d-1*(y4_3_3 - 1.0)*(y4_3_3 - 1.0) + 5.0d-1*(z4_3_3 - 1.0)*(z4_3_3 - 1.0) + 
	5.0d-1*(x5_3_3 - 1.0)*(x5_3_3 - 1.0) + 5.0d-1*(y5_3_3 - 1.0)*(y5_3_3 - 1.0) + 
	5.0d-1*(z5_3_3 - 1.0)*(z5_3_3 - 1.0) + 5.0d-1*(x6_3_3 - 1.0)*(x6_3_3 - 1.0) + 
	5.0d-1*(y6_3_3 - 1.0)*(y6_3_3 - 1.0) + 5.0d-1*(z6_3_3 - 1.0)*(z6_3_3 - 1.0) + 
	5.0d-1*(x7_3_3 - 1.0)*(x7_3_3 - 1.0) + 5.0d-1*(y7_3_3 - 1.0)*(y7_3_3 - 1.0) + 
	5.0d-1*(z7_3_3 - 1.0)*(z7_3_3 - 1.0) + 5.0d-1*(x8_3_3 - 1.0)*(x8_3_3 - 1.0) + 
	5.0d-1*(y8_3_3 - 1.0)*(y8_3_3 - 1.0) + 5.0d-1*(z8_3_3 - 1.0)*(z8_3_3 - 1.0) + 
	5.0d-1*(x9_3_3 - 1.0)*(x9_3_3 - 1.0) + 5.0d-1*(y9_3_3 - 1.0)*(y9_3_3 - 1.0) + 
	5.0d-1*(z9_3_3 - 1.0)*(z9_3_3 - 1.0) + 5.0d-1*(x1_4_3 - 1.0)*(x1_4_3 - 1.0) + 
	5.0d-1*(y1_4_3 - 1.0)*(y1_4_3 - 1.0) + 5.0d-1*(z1_4_3 - 1.0)*(z1_4_3 - 1.0) + 
	5.0d-1*(x2_4_3 - 1.0)*(x2_4_3 - 1.0) + 5.0d-1*(y2_4_3 - 1.0)*(y2_4_3 - 1.0) + 
	5.0d-1*(z2_4_3 - 1.0)*(z2_4_3 - 1.0) + 5.0d-1*(x3_4_3 - 1.0)*(x3_4_3 - 1.0) + 
	5.0d-1*(y3_4_3 - 1.0)*(y3_4_3 - 1.0) + 5.0d-1*(z3_4_3 - 1.0)*(z3_4_3 - 1.0) + 
	5.0d-1*(x4_4_3 - 1.0)*(x4_4_3 - 1.0) + 5.0d-1*(y4_4_3 - 1.0)*(y4_4_3 - 1.0) + 
	5.0d-1*(z4_4_3 - 1.0)*(z4_4_3 - 1.0) + 5.0d-1*(x5_4_3 - 1.0)*(x5_4_3 - 1.0) + 
	5.0d-1*(y5_4_3 - 1.0)*(y5_4_3 - 1.0) + 5.0d-1*(z5_4_3 - 1.0)*(z5_4_3 - 1.0) + 
	5.0d-1*(x6_4_3 - 1.0)*(x6_4_3 - 1.0) + 5.0d-1*(y6_4_3 - 1.0)*(y6_4_3 - 1.0) + 
	5.0d-1*(z6_4_3 - 1.0)*(z6_4_3 - 1.0) + 5.0d-1*(x7_4_3 - 1.0)*(x7_4_3 - 1.0) + 
	5.0d-1*(y7_4_3 - 1.0)*(y7_4_3 - 1.0) + 5.0d-1*(z7_4_3 - 1.0)*(z7_4_3 - 1.0) + 
	5.0d-1*(x8_4_3 - 1.0)*(x8_4_3 - 1.0) + 5.0d-1*(y8_4_3 - 1.0)*(y8_4_3 - 1.0) + 
	5.0d-1*(z8_4_3 - 1.0)*(z8_4_3 - 1.0) + 5.0d-1*(x9_4_3 - 1.0)*(x9_4_3 - 1.0) + 
	5.0d-1*(y9_4_3 - 1.0)*(y9_4_3 - 1.0) + 5.0d-1*(z9_4_3 - 1.0)*(z9_4_3 - 1.0) + 
	5.0d-1*(x1_5_3 - 1.0)*(x1_5_3 - 1.0) + 5.0d-1*(y1_5_3 - 1.0)*(y1_5_3 - 1.0) + 
	5.0d-1*(z1_5_3 - 1.0)*(z1_5_3 - 1.0) + 5.0d-1*(x2_5_3 - 1.0)*(x2_5_3 - 1.0) + 
	5.0d-1*(y2_5_3 - 1.0)*(y2_5_3 - 1.0) + 5.0d-1*(z2_5_3 - 1.0)*(z2_5_3 - 1.0) + 
	5.0d-1*(x3_5_3 - 1.0)*(x3_5_3 - 1.0) + 5.0d-1*(y3_5_3 - 1.0)*(y3_5_3 - 1.0) + 
	5.0d-1*(z3_5_3 - 1.0)*(z3_5_3 - 1.0) + 5.0d-1*(x4_5_3 - 1.0)*(x4_5_3 - 1.0) + 
	5.0d-1*(y4_5_3 - 1.0)*(y4_5_3 - 1.0) + 5.0d-1*(z4_5_3 - 1.0)*(z4_5_3 - 1.0) + 
	5.0d-1*(x5_5_3 - 1.0)*(x5_5_3 - 1.0) + 5.0d-1*(y5_5_3 - 1.0)*(y5_5_3 - 1.0) + 
	5.0d-1*(z5_5_3 - 1.0)*(z5_5_3 - 1.0) + 5.0d-1*(x6_5_3 - 1.0)*(x6_5_3 - 1.0) + 
	5.0d-1*(y6_5_3 - 1.0)*(y6_5_3 - 1.0) + 5.0d-1*(z6_5_3 - 1.0)*(z6_5_3 - 1.0) + 
	5.0d-1*(x7_5_3 - 1.0)*(x7_5_3 - 1.0) + 5.0d-1*(y7_5_3 - 1.0)*(y7_5_3 - 1.0) + 
	5.0d-1*(z7_5_3 - 1.0)*(z7_5_3 - 1.0) + 5.0d-1*(x8_5_3 - 1.0)*(x8_5_3 - 1.0) + 
	5.0d-1*(y8_5_3 - 1.0)*(y8_5_3 - 1.0) + 5.0d-1*(z8_5_3 - 1.0)*(z8_5_3 - 1.0) + 
	5.0d-1*(x9_5_3 - 1.0)*(x9_5_3 - 1.0) + 5.0d-1*(y9_5_3 - 1.0)*(y9_5_3 - 1.0) + 
	5.0d-1*(z9_5_3 - 1.0)*(z9_5_3 - 1.0) + 5.0d-1*(x1_6_3 - 1.0)*(x1_6_3 - 1.0) + 
	5.0d-1*(y1_6_3 - 1.0)*(y1_6_3 - 1.0) + 5.0d-1*(z1_6_3 - 1.0)*(z1_6_3 - 1.0) + 
	5.0d-1*(x2_6_3 - 1.0)*(x2_6_3 - 1.0) + 5.0d-1*(y2_6_3 - 1.0)*(y2_6_3 - 1.0) + 
	5.0d-1*(z2_6_3 - 1.0)*(z2_6_3 - 1.0) + 5.0d-1*(x3_6_3 - 1.0)*(x3_6_3 - 1.0) + 
	5.0d-1*(y3_6_3 - 1.0)*(y3_6_3 - 1.0) + 5.0d-1*(z3_6_3 - 1.0)*(z3_6_3 - 1.0) + 
	5.0d-1*(x4_6_3 - 1.0)*(x4_6_3 - 1.0) + 5.0d-1*(y4_6_3 - 1.0)*(y4_6_3 - 1.0) + 
	5.0d-1*(z4_6_3 - 1.0)*(z4_6_3 - 1.0) + 5.0d-1*(x5_6_3 - 1.0)*(x5_6_3 - 1.0) + 
	5.0d-1*(y5_6_3 - 1.0)*(y5_6_3 - 1.0) + 5.0d-1*(z5_6_3 - 1.0)*(z5_6_3 - 1.0) + 
	5.0d-1*(x6_6_3 - 1.0)*(x6_6_3 - 1.0) + 5.0d-1*(y6_6_3 - 1.0)*(y6_6_3 - 1.0) + 
	5.0d-1*(z6_6_3 - 1.0)*(z6_6_3 - 1.0) + 5.0d-1*(x7_6_3 - 1.0)*(x7_6_3 - 1.0) + 
	5.0d-1*(y7_6_3 - 1.0)*(y7_6_3 - 1.0) + 5.0d-1*(z7_6_3 - 1.0)*(z7_6_3 - 1.0) + 
	5.0d-1*(x8_6_3 - 1.0)*(x8_6_3 - 1.0) + 5.0d-1*(y8_6_3 - 1.0)*(y8_6_3 - 1.0) + 
	5.0d-1*(z8_6_3 - 1.0)*(z8_6_3 - 1.0) + 5.0d-1*(x9_6_3 - 1.0)*(x9_6_3 - 1.0) + 
	5.0d-1*(y9_6_3 - 1.0)*(y9_6_3 - 1.0) + 5.0d-1*(z9_6_3 - 1.0)*(z9_6_3 - 1.0) + 
	5.0d-1*(x1_7_3 - 1.0)*(x1_7_3 - 1.0) + 5.0d-1*(y1_7_3 - 1.0)*(y1_7_3 - 1.0) + 
	5.0d-1*(z1_7_3 - 1.0)*(z1_7_3 - 1.0) + 5.0d-1*(x2_7_3 - 1.0)*(x2_7_3 - 1.0) + 
	5.0d-1*(y2_7_3 - 1.0)*(y2_7_3 - 1.0) + 5.0d-1*(z2_7_3 - 1.0)*(z2_7_3 - 1.0) + 
	5.0d-1*(x3_7_3 - 1.0)*(x3_7_3 - 1.0) + 5.0d-1*(y3_7_3 - 1.0)*(y3_7_3 - 1.0) + 
	5.0d-1*(z3_7_3 - 1.0)*(z3_7_3 - 1.0) + 5.0d-1*(x4_7_3 - 1.0)*(x4_7_3 - 1.0) + 
	5.0d-1*(y4_7_3 - 1.0)*(y4_7_3 - 1.0) + 5.0d-1*(z4_7_3 - 1.0)*(z4_7_3 - 1.0) + 
	5.0d-1*(x5_7_3 - 1.0)*(x5_7_3 - 1.0) + 5.0d-1*(y5_7_3 - 1.0)*(y5_7_3 - 1.0) + 
	5.0d-1*(z5_7_3 - 1.0)*(z5_7_3 - 1.0) + 5.0d-1*(x6_7_3 - 1.0)*(x6_7_3 - 1.0) + 
	5.0d-1*(y6_7_3 - 1.0)*(y6_7_3 - 1.0) + 5.0d-1*(z6_7_3 - 1.0)*(z6_7_3 - 1.0) + 
	5.0d-1*(x7_7_3 - 1.0)*(x7_7_3 - 1.0) + 5.0d-1*(y7_7_3 - 1.0)*(y7_7_3 - 1.0) + 
	5.0d-1*(z7_7_3 - 1.0)*(z7_7_3 - 1.0) + 5.0d-1*(x8_7_3 - 1.0)*(x8_7_3 - 1.0) + 
	5.0d-1*(y8_7_3 - 1.0)*(y8_7_3 - 1.0) + 5.0d-1*(z8_7_3 - 1.0)*(z8_7_3 - 1.0) + 
	5.0d-1*(x9_7_3 - 1.0)*(x9_7_3 - 1.0) + 5.0d-1*(y9_7_3 - 1.0)*(y9_7_3 - 1.0) + 
	5.0d-1*(z9_7_3 - 1.0)*(z9_7_3 - 1.0) + 5.0d-1*(x1_8_3 - 1.0)*(x1_8_3 - 1.0) + 
	5.0d-1*(y1_8_3 - 1.0)*(y1_8_3 - 1.0) + 5.0d-1*(z1_8_3 - 1.0)*(z1_8_3 - 1.0) + 
	5.0d-1*(x2_8_3 - 1.0)*(x2_8_3 - 1.0) + 5.0d-1*(y2_8_3 - 1.0)*(y2_8_3 - 1.0) + 
	5.0d-1*(z2_8_3 - 1.0)*(z2_8_3 - 1.0) + 5.0d-1*(x3_8_3 - 1.0)*(x3_8_3 - 1.0) + 
	5.0d-1*(y3_8_3 - 1.0)*(y3_8_3 - 1.0) + 5.0d-1*(z3_8_3 - 1.0)*(z3_8_3 - 1.0) + 
	5.0d-1*(x4_8_3 - 1.0)*(x4_8_3 - 1.0) + 5.0d-1*(y4_8_3 - 1.0)*(y4_8_3 - 1.0) + 
	5.0d-1*(z4_8_3 - 1.0)*(z4_8_3 - 1.0) + 5.0d-1*(x5_8_3 - 1.0)*(x5_8_3 - 1.0) + 
	5.0d-1*(y5_8_3 - 1.0)*(y5_8_3 - 1.0) + 5.0d-1*(z5_8_3 - 1.0)*(z5_8_3 - 1.0) + 
	5.0d-1*(x6_8_3 - 1.0)*(x6_8_3 - 1.0) + 5.0d-1*(y6_8_3 - 1.0)*(y6_8_3 - 1.0) + 
	5.0d-1*(z6_8_3 - 1.0)*(z6_8_3 - 1.0) + 5.0d-1*(x7_8_3 - 1.0)*(x7_8_3 - 1.0) + 
	5.0d-1*(y7_8_3 - 1.0)*(y7_8_3 - 1.0) + 5.0d-1*(z7_8_3 - 1.0)*(z7_8_3 - 1.0) + 
	5.0d-1*(x8_8_3 - 1.0)*(x8_8_3 - 1.0) + 5.0d-1*(y8_8_3 - 1.0)*(y8_8_3 - 1.0) + 
	5.0d-1*(z8_8_3 - 1.0)*(z8_8_3 - 1.0) + 5.0d-1*(x9_8_3 - 1.0)*(x9_8_3 - 1.0) + 
	5.0d-1*(y9_8_3 - 1.0)*(y9_8_3 - 1.0) + 5.0d-1*(z9_8_3 - 1.0)*(z9_8_3 - 1.0) + 
	5.0d-1*(x1_9_3 - 1.0)*(x1_9_3 - 1.0) + 5.0d-1*(y1_9_3 - 1.0)*(y1_9_3 - 1.0) + 
	5.0d-1*(z1_9_3 - 1.0)*(z1_9_3 - 1.0) + 5.0d-1*(x2_9_3 - 1.0)*(x2_9_3 - 1.0) + 
	5.0d-1*(y2_9_3 - 1.0)*(y2_9_3 - 1.0) + 5.0d-1*(z2_9_3 - 1.0)*(z2_9_3 - 1.0) + 
	5.0d-1*(x3_9_3 - 1.0)*(x3_9_3 - 1.0) + 5.0d-1*(y3_9_3 - 1.0)*(y3_9_3 - 1.0) + 
	5.0d-1*(z3_9_3 - 1.0)*(z3_9_3 - 1.0) + 5.0d-1*(x4_9_3 - 1.0)*(x4_9_3 - 1.0) + 
	5.0d-1*(y4_9_3 - 1.0)*(y4_9_3 - 1.0) + 5.0d-1*(z4_9_3 - 1.0)*(z4_9_3 - 1.0) + 
	5.0d-1*(x5_9_3 - 1.0)*(x5_9_3 - 1.0) + 5.0d-1*(y5_9_3 - 1.0)*(y5_9_3 - 1.0) + 
	5.0d-1*(z5_9_3 - 1.0)*(z5_9_3 - 1.0) + 5.0d-1*(x6_9_3 - 1.0)*(x6_9_3 - 1.0) + 
	5.0d-1*(y6_9_3 - 1.0)*(y6_9_3 - 1.0) + 5.0d-1*(z6_9_3 - 1.0)*(z6_9_3 - 1.0) + 
	5.0d-1*(x7_9_3 - 1.0)*(x7_9_3 - 1.0) + 5.0d-1*(y7_9_3 - 1.0)*(y7_9_3 - 1.0) + 
	5.0d-1*(z7_9_3 - 1.0)*(z7_9_3 - 1.0) + 5.0d-1*(x8_9_3 - 1.0)*(x8_9_3 - 1.0) + 
	5.0d-1*(y8_9_3 - 1.0)*(y8_9_3 - 1.0) + 5.0d-1*(z8_9_3 - 1.0)*(z8_9_3 - 1.0) + 
	5.0d-1*(x9_9_3 - 1.0)*(x9_9_3 - 1.0) + 5.0d-1*(y9_9_3 - 1.0)*(y9_9_3 - 1.0) + 
	5.0d-1*(z9_9_3 - 1.0)*(z9_9_3 - 1.0) + 5.0d-1*(x1_1_4 - 1.0)*(x1_1_4 - 1.0) + 
	5.0d-1*(y1_1_4 - 1.0)*(y1_1_4 - 1.0) + 5.0d-1*(z1_1_4 - 1.0)*(z1_1_4 - 1.0) + 
	5.0d-1*(x2_1_4 - 1.0)*(x2_1_4 - 1.0) + 5.0d-1*(y2_1_4 - 1.0)*(y2_1_4 - 1.0) + 
	5.0d-1*(z2_1_4 - 1.0)*(z2_1_4 - 1.0) + 5.0d-1*(x3_1_4 - 1.0)*(x3_1_4 - 1.0) + 
	5.0d-1*(y3_1_4 - 1.0)*(y3_1_4 - 1.0) + 5.0d-1*(z3_1_4 - 1.0)*(z3_1_4 - 1.0) + 
	5.0d-1*(x4_1_4 - 1.0)*(x4_1_4 - 1.0) + 5.0d-1*(y4_1_4 - 1.0)*(y4_1_4 - 1.0) + 
	5.0d-1*(z4_1_4 - 1.0)*(z4_1_4 - 1.0) + 5.0d-1*(x5_1_4 - 1.0)*(x5_1_4 - 1.0) + 
	5.0d-1*(y5_1_4 - 1.0)*(y5_1_4 - 1.0) + 5.0d-1*(z5_1_4 - 1.0)*(z5_1_4 - 1.0) + 
	5.0d-1*(x6_1_4 - 1.0)*(x6_1_4 - 1.0) + 5.0d-1*(y6_1_4 - 1.0)*(y6_1_4 - 1.0) + 
	5.0d-1*(z6_1_4 - 1.0)*(z6_1_4 - 1.0) + 5.0d-1*(x7_1_4 - 1.0)*(x7_1_4 - 1.0) + 
	5.0d-1*(y7_1_4 - 1.0)*(y7_1_4 - 1.0) + 5.0d-1*(z7_1_4 - 1.0)*(z7_1_4 - 1.0) + 
	5.0d-1*(x8_1_4 - 1.0)*(x8_1_4 - 1.0) + 5.0d-1*(y8_1_4 - 1.0)*(y8_1_4 - 1.0) + 
	5.0d-1*(z8_1_4 - 1.0)*(z8_1_4 - 1.0) + 5.0d-1*(x9_1_4 - 1.0)*(x9_1_4 - 1.0) + 
	5.0d-1*(y9_1_4 - 1.0)*(y9_1_4 - 1.0) + 5.0d-1*(z9_1_4 - 1.0)*(z9_1_4 - 1.0) + 
	5.0d-1*(x1_2_4 - 1.0)*(x1_2_4 - 1.0) + 5.0d-1*(y1_2_4 - 1.0)*(y1_2_4 - 1.0) + 
	5.0d-1*(z1_2_4 - 1.0)*(z1_2_4 - 1.0) + 5.0d-1*(x2_2_4 - 1.0)*(x2_2_4 - 1.0) + 
	5.0d-1*(y2_2_4 - 1.0)*(y2_2_4 - 1.0) + 5.0d-1*(z2_2_4 - 1.0)*(z2_2_4 - 1.0) + 
	5.0d-1*(x3_2_4 - 1.0)*(x3_2_4 - 1.0) + 5.0d-1*(y3_2_4 - 1.0)*(y3_2_4 - 1.0) + 
	5.0d-1*(z3_2_4 - 1.0)*(z3_2_4 - 1.0) + 5.0d-1*(x4_2_4 - 1.0)*(x4_2_4 - 1.0) + 
	5.0d-1*(y4_2_4 - 1.0)*(y4_2_4 - 1.0) + 5.0d-1*(z4_2_4 - 1.0)*(z4_2_4 - 1.0) + 
	5.0d-1*(x5_2_4 - 1.0)*(x5_2_4 - 1.0) + 5.0d-1*(y5_2_4 - 1.0)*(y5_2_4 - 1.0) + 
	5.0d-1*(z5_2_4 - 1.0)*(z5_2_4 - 1.0) + 5.0d-1*(x6_2_4 - 1.0)*(x6_2_4 - 1.0) + 
	5.0d-1*(y6_2_4 - 1.0)*(y6_2_4 - 1.0) + 5.0d-1*(z6_2_4 - 1.0)*(z6_2_4 - 1.0) + 
	5.0d-1*(x7_2_4 - 1.0)*(x7_2_4 - 1.0) + 5.0d-1*(y7_2_4 - 1.0)*(y7_2_4 - 1.0) + 
	5.0d-1*(z7_2_4 - 1.0)*(z7_2_4 - 1.0) + 5.0d-1*(x8_2_4 - 1.0)*(x8_2_4 - 1.0) + 
	5.0d-1*(y8_2_4 - 1.0)*(y8_2_4 - 1.0) + 5.0d-1*(z8_2_4 - 1.0)*(z8_2_4 - 1.0) + 
	5.0d-1*(x9_2_4 - 1.0)*(x9_2_4 - 1.0) + 5.0d-1*(y9_2_4 - 1.0)*(y9_2_4 - 1.0) + 
	5.0d-1*(z9_2_4 - 1.0)*(z9_2_4 - 1.0) + 5.0d-1*(x1_3_4 - 1.0)*(x1_3_4 - 1.0) + 
	5.0d-1*(y1_3_4 - 1.0)*(y1_3_4 - 1.0) + 5.0d-1*(z1_3_4 - 1.0)*(z1_3_4 - 1.0) + 
	5.0d-1*(x2_3_4 - 1.0)*(x2_3_4 - 1.0) + 5.0d-1*(y2_3_4 - 1.0)*(y2_3_4 - 1.0) + 
	5.0d-1*(z2_3_4 - 1.0)*(z2_3_4 - 1.0) + 5.0d-1*(x3_3_4 - 1.0)*(x3_3_4 - 1.0) + 
	5.0d-1*(y3_3_4 - 1.0)*(y3_3_4 - 1.0) + 5.0d-1*(z3_3_4 - 1.0)*(z3_3_4 - 1.0) + 
	5.0d-1*(x4_3_4 - 1.0)*(x4_3_4 - 1.0) + 5.0d-1*(y4_3_4 - 1.0)*(y4_3_4 - 1.0) + 
	5.0d-1*(z4_3_4 - 1.0)*(z4_3_4 - 1.0) + 5.0d-1*(x5_3_4 - 1.0)*(x5_3_4 - 1.0) + 
	5.0d-1*(y5_3_4 - 1.0)*(y5_3_4 - 1.0) + 5.0d-1*(z5_3_4 - 1.0)*(z5_3_4 - 1.0) + 
	5.0d-1*(x6_3_4 - 1.0)*(x6_3_4 - 1.0) + 5.0d-1*(y6_3_4 - 1.0)*(y6_3_4 - 1.0) + 
	5.0d-1*(z6_3_4 - 1.0)*(z6_3_4 - 1.0) + 5.0d-1*(x7_3_4 - 1.0)*(x7_3_4 - 1.0) + 
	5.0d-1*(y7_3_4 - 1.0)*(y7_3_4 - 1.0) + 5.0d-1*(z7_3_4 - 1.0)*(z7_3_4 - 1.0) + 
	5.0d-1*(x8_3_4 - 1.0)*(x8_3_4 - 1.0) + 5.0d-1*(y8_3_4 - 1.0)*(y8_3_4 - 1.0) + 
	5.0d-1*(z8_3_4 - 1.0)*(z8_3_4 - 1.0) + 5.0d-1*(x9_3_4 - 1.0)*(x9_3_4 - 1.0) + 
	5.0d-1*(y9_3_4 - 1.0)*(y9_3_4 - 1.0) + 5.0d-1*(z9_3_4 - 1.0)*(z9_3_4 - 1.0) + 
	5.0d-1*(x1_4_4 - 1.0)*(x1_4_4 - 1.0) + 5.0d-1*(y1_4_4 - 1.0)*(y1_4_4 - 1.0) + 
	5.0d-1*(z1_4_4 - 1.0)*(z1_4_4 - 1.0) + 5.0d-1*(x2_4_4 - 1.0)*(x2_4_4 - 1.0) + 
	5.0d-1*(y2_4_4 - 1.0)*(y2_4_4 - 1.0) + 5.0d-1*(z2_4_4 - 1.0)*(z2_4_4 - 1.0) + 
	5.0d-1*(x3_4_4 - 1.0)*(x3_4_4 - 1.0) + 5.0d-1*(y3_4_4 - 1.0)*(y3_4_4 - 1.0) + 
	5.0d-1*(z3_4_4 - 1.0)*(z3_4_4 - 1.0) + 5.0d-1*(x4_4_4 - 1.0)*(x4_4_4 - 1.0) + 
	5.0d-1*(y4_4_4 - 1.0)*(y4_4_4 - 1.0) + 5.0d-1*(z4_4_4 - 1.0)*(z4_4_4 - 1.0) + 
	5.0d-1*(x5_4_4 - 1.0)*(x5_4_4 - 1.0) + 5.0d-1*(y5_4_4 - 1.0)*(y5_4_4 - 1.0) + 
	5.0d-1*(z5_4_4 - 1.0)*(z5_4_4 - 1.0) + 5.0d-1*(x6_4_4 - 1.0)*(x6_4_4 - 1.0) + 
	5.0d-1*(y6_4_4 - 1.0)*(y6_4_4 - 1.0) + 5.0d-1*(z6_4_4 - 1.0)*(z6_4_4 - 1.0) + 
	5.0d-1*(x7_4_4 - 1.0)*(x7_4_4 - 1.0) + 5.0d-1*(y7_4_4 - 1.0)*(y7_4_4 - 1.0) + 
	5.0d-1*(z7_4_4 - 1.0)*(z7_4_4 - 1.0) + 5.0d-1*(x8_4_4 - 1.0)*(x8_4_4 - 1.0) + 
	5.0d-1*(y8_4_4 - 1.0)*(y8_4_4 - 1.0) + 5.0d-1*(z8_4_4 - 1.0)*(z8_4_4 - 1.0) + 
	5.0d-1*(x9_4_4 - 1.0)*(x9_4_4 - 1.0) + 5.0d-1*(y9_4_4 - 1.0)*(y9_4_4 - 1.0) + 
	5.0d-1*(z9_4_4 - 1.0)*(z9_4_4 - 1.0) + 5.0d-1*(x1_5_4 - 1.0)*(x1_5_4 - 1.0) + 
	5.0d-1*(y1_5_4 - 1.0)*(y1_5_4 - 1.0) + 5.0d-1*(z1_5_4 - 1.0)*(z1_5_4 - 1.0) + 
	5.0d-1*(x2_5_4 - 1.0)*(x2_5_4 - 1.0) + 5.0d-1*(y2_5_4 - 1.0)*(y2_5_4 - 1.0) + 
	5.0d-1*(z2_5_4 - 1.0)*(z2_5_4 - 1.0) + 5.0d-1*(x3_5_4 - 1.0)*(x3_5_4 - 1.0) + 
	5.0d-1*(y3_5_4 - 1.0)*(y3_5_4 - 1.0) + 5.0d-1*(z3_5_4 - 1.0)*(z3_5_4 - 1.0) + 
	5.0d-1*(x4_5_4 - 1.0)*(x4_5_4 - 1.0) + 5.0d-1*(y4_5_4 - 1.0)*(y4_5_4 - 1.0) + 
	5.0d-1*(z4_5_4 - 1.0)*(z4_5_4 - 1.0) + 5.0d-1*(x5_5_4 - 1.0)*(x5_5_4 - 1.0) + 
	5.0d-1*(y5_5_4 - 1.0)*(y5_5_4 - 1.0) + 5.0d-1*(z5_5_4 - 1.0)*(z5_5_4 - 1.0) + 
	5.0d-1*(x6_5_4 - 1.0)*(x6_5_4 - 1.0) + 5.0d-1*(y6_5_4 - 1.0)*(y6_5_4 - 1.0) + 
	5.0d-1*(z6_5_4 - 1.0)*(z6_5_4 - 1.0) + 5.0d-1*(x7_5_4 - 1.0)*(x7_5_4 - 1.0) + 
	5.0d-1*(y7_5_4 - 1.0)*(y7_5_4 - 1.0) + 5.0d-1*(z7_5_4 - 1.0)*(z7_5_4 - 1.0) + 
	5.0d-1*(x8_5_4 - 1.0)*(x8_5_4 - 1.0) + 5.0d-1*(y8_5_4 - 1.0)*(y8_5_4 - 1.0) + 
	5.0d-1*(z8_5_4 - 1.0)*(z8_5_4 - 1.0) + 5.0d-1*(x9_5_4 - 1.0)*(x9_5_4 - 1.0) + 
	5.0d-1*(y9_5_4 - 1.0)*(y9_5_4 - 1.0) + 5.0d-1*(z9_5_4 - 1.0)*(z9_5_4 - 1.0) + 
	5.0d-1*(x1_6_4 - 1.0)*(x1_6_4 - 1.0) + 5.0d-1*(y1_6_4 - 1.0)*(y1_6_4 - 1.0) + 
	5.0d-1*(z1_6_4 - 1.0)*(z1_6_4 - 1.0) + 5.0d-1*(x2_6_4 - 1.0)*(x2_6_4 - 1.0) + 
	5.0d-1*(y2_6_4 - 1.0)*(y2_6_4 - 1.0) + 5.0d-1*(z2_6_4 - 1.0)*(z2_6_4 - 1.0) + 
	5.0d-1*(x3_6_4 - 1.0)*(x3_6_4 - 1.0) + 5.0d-1*(y3_6_4 - 1.0)*(y3_6_4 - 1.0) + 
	5.0d-1*(z3_6_4 - 1.0)*(z3_6_4 - 1.0) + 5.0d-1*(x4_6_4 - 1.0)*(x4_6_4 - 1.0) + 
	5.0d-1*(y4_6_4 - 1.0)*(y4_6_4 - 1.0) + 5.0d-1*(z4_6_4 - 1.0)*(z4_6_4 - 1.0) + 
	5.0d-1*(x5_6_4 - 1.0)*(x5_6_4 - 1.0) + 5.0d-1*(y5_6_4 - 1.0)*(y5_6_4 - 1.0) + 
	5.0d-1*(z5_6_4 - 1.0)*(z5_6_4 - 1.0) + 5.0d-1*(x6_6_4 - 1.0)*(x6_6_4 - 1.0) + 
	5.0d-1*(y6_6_4 - 1.0)*(y6_6_4 - 1.0) + 5.0d-1*(z6_6_4 - 1.0)*(z6_6_4 - 1.0) + 
	5.0d-1*(x7_6_4 - 1.0)*(x7_6_4 - 1.0) + 5.0d-1*(y7_6_4 - 1.0)*(y7_6_4 - 1.0) + 
	5.0d-1*(z7_6_4 - 1.0)*(z7_6_4 - 1.0) + 5.0d-1*(x8_6_4 - 1.0)*(x8_6_4 - 1.0) + 
	5.0d-1*(y8_6_4 - 1.0)*(y8_6_4 - 1.0) + 5.0d-1*(z8_6_4 - 1.0)*(z8_6_4 - 1.0) + 
	5.0d-1*(x9_6_4 - 1.0)*(x9_6_4 - 1.0) + 5.0d-1*(y9_6_4 - 1.0)*(y9_6_4 - 1.0) + 
	5.0d-1*(z9_6_4 - 1.0)*(z9_6_4 - 1.0) + 5.0d-1*(x1_7_4 - 1.0)*(x1_7_4 - 1.0) + 
	5.0d-1*(y1_7_4 - 1.0)*(y1_7_4 - 1.0) + 5.0d-1*(z1_7_4 - 1.0)*(z1_7_4 - 1.0) + 
	5.0d-1*(x2_7_4 - 1.0)*(x2_7_4 - 1.0) + 5.0d-1*(y2_7_4 - 1.0)*(y2_7_4 - 1.0) + 
	5.0d-1*(z2_7_4 - 1.0)*(z2_7_4 - 1.0) + 5.0d-1*(x3_7_4 - 1.0)*(x3_7_4 - 1.0) + 
	5.0d-1*(y3_7_4 - 1.0)*(y3_7_4 - 1.0) + 5.0d-1*(z3_7_4 - 1.0)*(z3_7_4 - 1.0) + 
	5.0d-1*(x4_7_4 - 1.0)*(x4_7_4 - 1.0) + 5.0d-1*(y4_7_4 - 1.0)*(y4_7_4 - 1.0) + 
	5.0d-1*(z4_7_4 - 1.0)*(z4_7_4 - 1.0) + 5.0d-1*(x5_7_4 - 1.0)*(x5_7_4 - 1.0) + 
	5.0d-1*(y5_7_4 - 1.0)*(y5_7_4 - 1.0) + 5.0d-1*(z5_7_4 - 1.0)*(z5_7_4 - 1.0) + 
	5.0d-1*(x6_7_4 - 1.0)*(x6_7_4 - 1.0) + 5.0d-1*(y6_7_4 - 1.0)*(y6_7_4 - 1.0) + 
	5.0d-1*(z6_7_4 - 1.0)*(z6_7_4 - 1.0) + 5.0d-1*(x7_7_4 - 1.0)*(x7_7_4 - 1.0) + 
	5.0d-1*(y7_7_4 - 1.0)*(y7_7_4 - 1.0) + 5.0d-1*(z7_7_4 - 1.0)*(z7_7_4 - 1.0) + 
	5.0d-1*(x8_7_4 - 1.0)*(x8_7_4 - 1.0) + 5.0d-1*(y8_7_4 - 1.0)*(y8_7_4 - 1.0) + 
	5.0d-1*(z8_7_4 - 1.0)*(z8_7_4 - 1.0) + 5.0d-1*(x9_7_4 - 1.0)*(x9_7_4 - 1.0) + 
	5.0d-1*(y9_7_4 - 1.0)*(y9_7_4 - 1.0) + 5.0d-1*(z9_7_4 - 1.0)*(z9_7_4 - 1.0) + 
	5.0d-1*(x1_8_4 - 1.0)*(x1_8_4 - 1.0) + 5.0d-1*(y1_8_4 - 1.0)*(y1_8_4 - 1.0) + 
	5.0d-1*(z1_8_4 - 1.0)*(z1_8_4 - 1.0) + 5.0d-1*(x2_8_4 - 1.0)*(x2_8_4 - 1.0) + 
	5.0d-1*(y2_8_4 - 1.0)*(y2_8_4 - 1.0) + 5.0d-1*(z2_8_4 - 1.0)*(z2_8_4 - 1.0) + 
	5.0d-1*(x3_8_4 - 1.0)*(x3_8_4 - 1.0) + 5.0d-1*(y3_8_4 - 1.0)*(y3_8_4 - 1.0) + 
	5.0d-1*(z3_8_4 - 1.0)*(z3_8_4 - 1.0) + 5.0d-1*(x4_8_4 - 1.0)*(x4_8_4 - 1.0) + 
	5.0d-1*(y4_8_4 - 1.0)*(y4_8_4 - 1.0) + 5.0d-1*(z4_8_4 - 1.0)*(z4_8_4 - 1.0) + 
	5.0d-1*(x5_8_4 - 1.0)*(x5_8_4 - 1.0) + 5.0d-1*(y5_8_4 - 1.0)*(y5_8_4 - 1.0) + 
	5.0d-1*(z5_8_4 - 1.0)*(z5_8_4 - 1.0) + 5.0d-1*(x6_8_4 - 1.0)*(x6_8_4 - 1.0) + 
	5.0d-1*(y6_8_4 - 1.0)*(y6_8_4 - 1.0) + 5.0d-1*(z6_8_4 - 1.0)*(z6_8_4 - 1.0) + 
	5.0d-1*(x7_8_4 - 1.0)*(x7_8_4 - 1.0) + 5.0d-1*(y7_8_4 - 1.0)*(y7_8_4 - 1.0) + 
	5.0d-1*(z7_8_4 - 1.0)*(z7_8_4 - 1.0) + 5.0d-1*(x8_8_4 - 1.0)*(x8_8_4 - 1.0) + 
	5.0d-1*(y8_8_4 - 1.0)*(y8_8_4 - 1.0) + 5.0d-1*(z8_8_4 - 1.0)*(z8_8_4 - 1.0) + 
	5.0d-1*(x9_8_4 - 1.0)*(x9_8_4 - 1.0) + 5.0d-1*(y9_8_4 - 1.0)*(y9_8_4 - 1.0) + 
	5.0d-1*(z9_8_4 - 1.0)*(z9_8_4 - 1.0) + 5.0d-1*(x1_9_4 - 1.0)*(x1_9_4 - 1.0) + 
	5.0d-1*(y1_9_4 - 1.0)*(y1_9_4 - 1.0) + 5.0d-1*(z1_9_4 - 1.0)*(z1_9_4 - 1.0) + 
	5.0d-1*(x2_9_4 - 1.0)*(x2_9_4 - 1.0) + 5.0d-1*(y2_9_4 - 1.0)*(y2_9_4 - 1.0) + 
	5.0d-1*(z2_9_4 - 1.0)*(z2_9_4 - 1.0) + 5.0d-1*(x3_9_4 - 1.0)*(x3_9_4 - 1.0) + 
	5.0d-1*(y3_9_4 - 1.0)*(y3_9_4 - 1.0) + 5.0d-1*(z3_9_4 - 1.0)*(z3_9_4 - 1.0) + 
	5.0d-1*(x4_9_4 - 1.0)*(x4_9_4 - 1.0) + 5.0d-1*(y4_9_4 - 1.0)*(y4_9_4 - 1.0) + 
	5.0d-1*(z4_9_4 - 1.0)*(z4_9_4 - 1.0) + 5.0d-1*(x5_9_4 - 1.0)*(x5_9_4 - 1.0) + 
	5.0d-1*(y5_9_4 - 1.0)*(y5_9_4 - 1.0) + 5.0d-1*(z5_9_4 - 1.0)*(z5_9_4 - 1.0) + 
	5.0d-1*(x6_9_4 - 1.0)*(x6_9_4 - 1.0) + 5.0d-1*(y6_9_4 - 1.0)*(y6_9_4 - 1.0) + 
	5.0d-1*(z6_9_4 - 1.0)*(z6_9_4 - 1.0) + 5.0d-1*(x7_9_4 - 1.0)*(x7_9_4 - 1.0) + 
	5.0d-1*(y7_9_4 - 1.0)*(y7_9_4 - 1.0) + 5.0d-1*(z7_9_4 - 1.0)*(z7_9_4 - 1.0) + 
	5.0d-1*(x8_9_4 - 1.0)*(x8_9_4 - 1.0) + 5.0d-1*(y8_9_4 - 1.0)*(y8_9_4 - 1.0) + 
	5.0d-1*(z8_9_4 - 1.0)*(z8_9_4 - 1.0) + 5.0d-1*(x9_9_4 - 1.0)*(x9_9_4 - 1.0) + 
	5.0d-1*(y9_9_4 - 1.0)*(y9_9_4 - 1.0) + 5.0d-1*(z9_9_4 - 1.0)*(z9_9_4 - 1.0) + 
	5.0d-1*(x1_1_5 - 1.0)*(x1_1_5 - 1.0) + 5.0d-1*(y1_1_5 - 1.0)*(y1_1_5 - 1.0) + 
	5.0d-1*(z1_1_5 - 1.0)*(z1_1_5 - 1.0) + 5.0d-1*(x2_1_5 - 1.0)*(x2_1_5 - 1.0) + 
	5.0d-1*(y2_1_5 - 1.0)*(y2_1_5 - 1.0) + 5.0d-1*(z2_1_5 - 1.0)*(z2_1_5 - 1.0) + 
	5.0d-1*(x3_1_5 - 1.0)*(x3_1_5 - 1.0) + 5.0d-1*(y3_1_5 - 1.0)*(y3_1_5 - 1.0) + 
	5.0d-1*(z3_1_5 - 1.0)*(z3_1_5 - 1.0) + 5.0d-1*(x4_1_5 - 1.0)*(x4_1_5 - 1.0) + 
	5.0d-1*(y4_1_5 - 1.0)*(y4_1_5 - 1.0) + 5.0d-1*(z4_1_5 - 1.0)*(z4_1_5 - 1.0) + 
	5.0d-1*(x5_1_5 - 1.0)*(x5_1_5 - 1.0) + 5.0d-1*(y5_1_5 - 1.0)*(y5_1_5 - 1.0) + 
	5.0d-1*(z5_1_5 - 1.0)*(z5_1_5 - 1.0) + 5.0d-1*(x6_1_5 - 1.0)*(x6_1_5 - 1.0) + 
	5.0d-1*(y6_1_5 - 1.0)*(y6_1_5 - 1.0) + 5.0d-1*(z6_1_5 - 1.0)*(z6_1_5 - 1.0) + 
	5.0d-1*(x7_1_5 - 1.0)*(x7_1_5 - 1.0) + 5.0d-1*(y7_1_5 - 1.0)*(y7_1_5 - 1.0) + 
	5.0d-1*(z7_1_5 - 1.0)*(z7_1_5 - 1.0) + 5.0d-1*(x8_1_5 - 1.0)*(x8_1_5 - 1.0) + 
	5.0d-1*(y8_1_5 - 1.0)*(y8_1_5 - 1.0) + 5.0d-1*(z8_1_5 - 1.0)*(z8_1_5 - 1.0) + 
	5.0d-1*(x9_1_5 - 1.0)*(x9_1_5 - 1.0) + 5.0d-1*(y9_1_5 - 1.0)*(y9_1_5 - 1.0) + 
	5.0d-1*(z9_1_5 - 1.0)*(z9_1_5 - 1.0) + 5.0d-1*(x1_2_5 - 1.0)*(x1_2_5 - 1.0) + 
	5.0d-1*(y1_2_5 - 1.0)*(y1_2_5 - 1.0) + 5.0d-1*(z1_2_5 - 1.0)*(z1_2_5 - 1.0) + 
	5.0d-1*(x2_2_5 - 1.0)*(x2_2_5 - 1.0) + 5.0d-1*(y2_2_5 - 1.0)*(y2_2_5 - 1.0) + 
	5.0d-1*(z2_2_5 - 1.0)*(z2_2_5 - 1.0) + 5.0d-1*(x3_2_5 - 1.0)*(x3_2_5 - 1.0) + 
	5.0d-1*(y3_2_5 - 1.0)*(y3_2_5 - 1.0) + 5.0d-1*(z3_2_5 - 1.0)*(z3_2_5 - 1.0) + 
	5.0d-1*(x4_2_5 - 1.0)*(x4_2_5 - 1.0) + 5.0d-1*(y4_2_5 - 1.0)*(y4_2_5 - 1.0) + 
	5.0d-1*(z4_2_5 - 1.0)*(z4_2_5 - 1.0) + 5.0d-1*(x5_2_5 - 1.0)*(x5_2_5 - 1.0) + 
	5.0d-1*(y5_2_5 - 1.0)*(y5_2_5 - 1.0) + 5.0d-1*(z5_2_5 - 1.0)*(z5_2_5 - 1.0) + 
	5.0d-1*(x6_2_5 - 1.0)*(x6_2_5 - 1.0) + 5.0d-1*(y6_2_5 - 1.0)*(y6_2_5 - 1.0) + 
	5.0d-1*(z6_2_5 - 1.0)*(z6_2_5 - 1.0) + 5.0d-1*(x7_2_5 - 1.0)*(x7_2_5 - 1.0) + 
	5.0d-1*(y7_2_5 - 1.0)*(y7_2_5 - 1.0) + 5.0d-1*(z7_2_5 - 1.0)*(z7_2_5 - 1.0) + 
	5.0d-1*(x8_2_5 - 1.0)*(x8_2_5 - 1.0) + 5.0d-1*(y8_2_5 - 1.0)*(y8_2_5 - 1.0) + 
	5.0d-1*(z8_2_5 - 1.0)*(z8_2_5 - 1.0) + 5.0d-1*(x9_2_5 - 1.0)*(x9_2_5 - 1.0) + 
	5.0d-1*(y9_2_5 - 1.0)*(y9_2_5 - 1.0) + 5.0d-1*(z9_2_5 - 1.0)*(z9_2_5 - 1.0) + 
	5.0d-1*(x1_3_5 - 1.0)*(x1_3_5 - 1.0) + 5.0d-1*(y1_3_5 - 1.0)*(y1_3_5 - 1.0) + 
	5.0d-1*(z1_3_5 - 1.0)*(z1_3_5 - 1.0) + 5.0d-1*(x2_3_5 - 1.0)*(x2_3_5 - 1.0) + 
	5.0d-1*(y2_3_5 - 1.0)*(y2_3_5 - 1.0) + 5.0d-1*(z2_3_5 - 1.0)*(z2_3_5 - 1.0) + 
	5.0d-1*(x3_3_5 - 1.0)*(x3_3_5 - 1.0) + 5.0d-1*(y3_3_5 - 1.0)*(y3_3_5 - 1.0) + 
	5.0d-1*(z3_3_5 - 1.0)*(z3_3_5 - 1.0) + 5.0d-1*(x4_3_5 - 1.0)*(x4_3_5 - 1.0) + 
	5.0d-1*(y4_3_5 - 1.0)*(y4_3_5 - 1.0) + 5.0d-1*(z4_3_5 - 1.0)*(z4_3_5 - 1.0) + 
	5.0d-1*(x5_3_5 - 1.0)*(x5_3_5 - 1.0) + 5.0d-1*(y5_3_5 - 1.0)*(y5_3_5 - 1.0) + 
	5.0d-1*(z5_3_5 - 1.0)*(z5_3_5 - 1.0) + 5.0d-1*(x6_3_5 - 1.0)*(x6_3_5 - 1.0) + 
	5.0d-1*(y6_3_5 - 1.0)*(y6_3_5 - 1.0) + 5.0d-1*(z6_3_5 - 1.0)*(z6_3_5 - 1.0) + 
	5.0d-1*(x7_3_5 - 1.0)*(x7_3_5 - 1.0) + 5.0d-1*(y7_3_5 - 1.0)*(y7_3_5 - 1.0) + 
	5.0d-1*(z7_3_5 - 1.0)*(z7_3_5 - 1.0) + 5.0d-1*(x8_3_5 - 1.0)*(x8_3_5 - 1.0) + 
	5.0d-1*(y8_3_5 - 1.0)*(y8_3_5 - 1.0) + 5.0d-1*(z8_3_5 - 1.0)*(z8_3_5 - 1.0) + 
	5.0d-1*(x9_3_5 - 1.0)*(x9_3_5 - 1.0) + 5.0d-1*(y9_3_5 - 1.0)*(y9_3_5 - 1.0) + 
	5.0d-1*(z9_3_5 - 1.0)*(z9_3_5 - 1.0) + 5.0d-1*(x1_4_5 - 1.0)*(x1_4_5 - 1.0) + 
	5.0d-1*(y1_4_5 - 1.0)*(y1_4_5 - 1.0) + 5.0d-1*(z1_4_5 - 1.0)*(z1_4_5 - 1.0) + 
	5.0d-1*(x2_4_5 - 1.0)*(x2_4_5 - 1.0) + 5.0d-1*(y2_4_5 - 1.0)*(y2_4_5 - 1.0) + 
	5.0d-1*(z2_4_5 - 1.0)*(z2_4_5 - 1.0) + 5.0d-1*(x3_4_5 - 1.0)*(x3_4_5 - 1.0) + 
	5.0d-1*(y3_4_5 - 1.0)*(y3_4_5 - 1.0) + 5.0d-1*(z3_4_5 - 1.0)*(z3_4_5 - 1.0) + 
	5.0d-1*(x4_4_5 - 1.0)*(x4_4_5 - 1.0) + 5.0d-1*(y4_4_5 - 1.0)*(y4_4_5 - 1.0) + 
	5.0d-1*(z4_4_5 - 1.0)*(z4_4_5 - 1.0) + 5.0d-1*(x5_4_5 - 1.0)*(x5_4_5 - 1.0) + 
	5.0d-1*(y5_4_5 - 1.0)*(y5_4_5 - 1.0) + 5.0d-1*(z5_4_5 - 1.0)*(z5_4_5 - 1.0) + 
	5.0d-1*(x6_4_5 - 1.0)*(x6_4_5 - 1.0) + 5.0d-1*(y6_4_5 - 1.0)*(y6_4_5 - 1.0) + 
	5.0d-1*(z6_4_5 - 1.0)*(z6_4_5 - 1.0) + 5.0d-1*(x7_4_5 - 1.0)*(x7_4_5 - 1.0) + 
	5.0d-1*(y7_4_5 - 1.0)*(y7_4_5 - 1.0) + 5.0d-1*(z7_4_5 - 1.0)*(z7_4_5 - 1.0) + 
	5.0d-1*(x8_4_5 - 1.0)*(x8_4_5 - 1.0) + 5.0d-1*(y8_4_5 - 1.0)*(y8_4_5 - 1.0) + 
	5.0d-1*(z8_4_5 - 1.0)*(z8_4_5 - 1.0) + 5.0d-1*(x9_4_5 - 1.0)*(x9_4_5 - 1.0) + 
	5.0d-1*(y9_4_5 - 1.0)*(y9_4_5 - 1.0) + 5.0d-1*(z9_4_5 - 1.0)*(z9_4_5 - 1.0) + 
	5.0d-1*(x1_5_5 - 1.0)*(x1_5_5 - 1.0) + 5.0d-1*(y1_5_5 - 1.0)*(y1_5_5 - 1.0) + 
	5.0d-1*(z1_5_5 - 1.0)*(z1_5_5 - 1.0) + 5.0d-1*(x2_5_5 - 1.0)*(x2_5_5 - 1.0) + 
	5.0d-1*(y2_5_5 - 1.0)*(y2_5_5 - 1.0) + 5.0d-1*(z2_5_5 - 1.0)*(z2_5_5 - 1.0) + 
	5.0d-1*(x3_5_5 - 1.0)*(x3_5_5 - 1.0) + 5.0d-1*(y3_5_5 - 1.0)*(y3_5_5 - 1.0) + 
	5.0d-1*(z3_5_5 - 1.0)*(z3_5_5 - 1.0) + 5.0d-1*(x4_5_5 - 1.0)*(x4_5_5 - 1.0) + 
	5.0d-1*(y4_5_5 - 1.0)*(y4_5_5 - 1.0) + 5.0d-1*(z4_5_5 - 1.0)*(z4_5_5 - 1.0) + 
	5.0d-1*(x5_5_5 - 1.0)*(x5_5_5 - 1.0) + 5.0d-1*(y5_5_5 - 1.0)*(y5_5_5 - 1.0) + 
	5.0d-1*(z5_5_5 - 1.0)*(z5_5_5 - 1.0) + 5.0d-1*(x6_5_5 - 1.0)*(x6_5_5 - 1.0) + 
	5.0d-1*(y6_5_5 - 1.0)*(y6_5_5 - 1.0) + 5.0d-1*(z6_5_5 - 1.0)*(z6_5_5 - 1.0) + 
	5.0d-1*(x7_5_5 - 1.0)*(x7_5_5 - 1.0) + 5.0d-1*(y7_5_5 - 1.0)*(y7_5_5 - 1.0) + 
	5.0d-1*(z7_5_5 - 1.0)*(z7_5_5 - 1.0) + 5.0d-1*(x8_5_5 - 1.0)*(x8_5_5 - 1.0) + 
	5.0d-1*(y8_5_5 - 1.0)*(y8_5_5 - 1.0) + 5.0d-1*(z8_5_5 - 1.0)*(z8_5_5 - 1.0) + 
	5.0d-1*(x9_5_5 - 1.0)*(x9_5_5 - 1.0) + 5.0d-1*(y9_5_5 - 1.0)*(y9_5_5 - 1.0) + 
	5.0d-1*(z9_5_5 - 1.0)*(z9_5_5 - 1.0) + 5.0d-1*(x1_6_5 - 1.0)*(x1_6_5 - 1.0) + 
	5.0d-1*(y1_6_5 - 1.0)*(y1_6_5 - 1.0) + 5.0d-1*(z1_6_5 - 1.0)*(z1_6_5 - 1.0) + 
	5.0d-1*(x2_6_5 - 1.0)*(x2_6_5 - 1.0) + 5.0d-1*(y2_6_5 - 1.0)*(y2_6_5 - 1.0) + 
	5.0d-1*(z2_6_5 - 1.0)*(z2_6_5 - 1.0) + 5.0d-1*(x3_6_5 - 1.0)*(x3_6_5 - 1.0) + 
	5.0d-1*(y3_6_5 - 1.0)*(y3_6_5 - 1.0) + 5.0d-1*(z3_6_5 - 1.0)*(z3_6_5 - 1.0) + 
	5.0d-1*(x4_6_5 - 1.0)*(x4_6_5 - 1.0) + 5.0d-1*(y4_6_5 - 1.0)*(y4_6_5 - 1.0) + 
	5.0d-1*(z4_6_5 - 1.0)*(z4_6_5 - 1.0) + 5.0d-1*(x5_6_5 - 1.0)*(x5_6_5 - 1.0) + 
	5.0d-1*(y5_6_5 - 1.0)*(y5_6_5 - 1.0) + 5.0d-1*(z5_6_5 - 1.0)*(z5_6_5 - 1.0) + 
	5.0d-1*(x6_6_5 - 1.0)*(x6_6_5 - 1.0) + 5.0d-1*(y6_6_5 - 1.0)*(y6_6_5 - 1.0) + 
	5.0d-1*(z6_6_5 - 1.0)*(z6_6_5 - 1.0) + 5.0d-1*(x7_6_5 - 1.0)*(x7_6_5 - 1.0) + 
	5.0d-1*(y7_6_5 - 1.0)*(y7_6_5 - 1.0) + 5.0d-1*(z7_6_5 - 1.0)*(z7_6_5 - 1.0) + 
	5.0d-1*(x8_6_5 - 1.0)*(x8_6_5 - 1.0) + 5.0d-1*(y8_6_5 - 1.0)*(y8_6_5 - 1.0) + 
	5.0d-1*(z8_6_5 - 1.0)*(z8_6_5 - 1.0) + 5.0d-1*(x9_6_5 - 1.0)*(x9_6_5 - 1.0) + 
	5.0d-1*(y9_6_5 - 1.0)*(y9_6_5 - 1.0) + 5.0d-1*(z9_6_5 - 1.0)*(z9_6_5 - 1.0) + 
	5.0d-1*(x1_7_5 - 1.0)*(x1_7_5 - 1.0) + 5.0d-1*(y1_7_5 - 1.0)*(y1_7_5 - 1.0) + 
	5.0d-1*(z1_7_5 - 1.0)*(z1_7_5 - 1.0) + 5.0d-1*(x2_7_5 - 1.0)*(x2_7_5 - 1.0) + 
	5.0d-1*(y2_7_5 - 1.0)*(y2_7_5 - 1.0) + 5.0d-1*(z2_7_5 - 1.0)*(z2_7_5 - 1.0) + 
	5.0d-1*(x3_7_5 - 1.0)*(x3_7_5 - 1.0) + 5.0d-1*(y3_7_5 - 1.0)*(y3_7_5 - 1.0) + 
	5.0d-1*(z3_7_5 - 1.0)*(z3_7_5 - 1.0) + 5.0d-1*(x4_7_5 - 1.0)*(x4_7_5 - 1.0) + 
	5.0d-1*(y4_7_5 - 1.0)*(y4_7_5 - 1.0) + 5.0d-1*(z4_7_5 - 1.0)*(z4_7_5 - 1.0) + 
	5.0d-1*(x5_7_5 - 1.0)*(x5_7_5 - 1.0) + 5.0d-1*(y5_7_5 - 1.0)*(y5_7_5 - 1.0) + 
	5.0d-1*(z5_7_5 - 1.0)*(z5_7_5 - 1.0) + 5.0d-1*(x6_7_5 - 1.0)*(x6_7_5 - 1.0) + 
	5.0d-1*(y6_7_5 - 1.0)*(y6_7_5 - 1.0) + 5.0d-1*(z6_7_5 - 1.0)*(z6_7_5 - 1.0) + 
	5.0d-1*(x7_7_5 - 1.0)*(x7_7_5 - 1.0) + 5.0d-1*(y7_7_5 - 1.0)*(y7_7_5 - 1.0) + 
	5.0d-1*(z7_7_5 - 1.0)*(z7_7_5 - 1.0) + 5.0d-1*(x8_7_5 - 1.0)*(x8_7_5 - 1.0) + 
	5.0d-1*(y8_7_5 - 1.0)*(y8_7_5 - 1.0) + 5.0d-1*(z8_7_5 - 1.0)*(z8_7_5 - 1.0) + 
	5.0d-1*(x9_7_5 - 1.0)*(x9_7_5 - 1.0) + 5.0d-1*(y9_7_5 - 1.0)*(y9_7_5 - 1.0) + 
	5.0d-1*(z9_7_5 - 1.0)*(z9_7_5 - 1.0) + 5.0d-1*(x1_8_5 - 1.0)*(x1_8_5 - 1.0) + 
	5.0d-1*(y1_8_5 - 1.0)*(y1_8_5 - 1.0) + 5.0d-1*(z1_8_5 - 1.0)*(z1_8_5 - 1.0) + 
	5.0d-1*(x2_8_5 - 1.0)*(x2_8_5 - 1.0) + 5.0d-1*(y2_8_5 - 1.0)*(y2_8_5 - 1.0) + 
	5.0d-1*(z2_8_5 - 1.0)*(z2_8_5 - 1.0) + 5.0d-1*(x3_8_5 - 1.0)*(x3_8_5 - 1.0) + 
	5.0d-1*(y3_8_5 - 1.0)*(y3_8_5 - 1.0) + 5.0d-1*(z3_8_5 - 1.0)*(z3_8_5 - 1.0) + 
	5.0d-1*(x4_8_5 - 1.0)*(x4_8_5 - 1.0) + 5.0d-1*(y4_8_5 - 1.0)*(y4_8_5 - 1.0) + 
	5.0d-1*(z4_8_5 - 1.0)*(z4_8_5 - 1.0) + 5.0d-1*(x5_8_5 - 1.0)*(x5_8_5 - 1.0) + 
	5.0d-1*(y5_8_5 - 1.0)*(y5_8_5 - 1.0) + 5.0d-1*(z5_8_5 - 1.0)*(z5_8_5 - 1.0) + 
	5.0d-1*(x6_8_5 - 1.0)*(x6_8_5 - 1.0) + 5.0d-1*(y6_8_5 - 1.0)*(y6_8_5 - 1.0) + 
	5.0d-1*(z6_8_5 - 1.0)*(z6_8_5 - 1.0) + 5.0d-1*(x7_8_5 - 1.0)*(x7_8_5 - 1.0) + 
	5.0d-1*(y7_8_5 - 1.0)*(y7_8_5 - 1.0) + 5.0d-1*(z7_8_5 - 1.0)*(z7_8_5 - 1.0) + 
	5.0d-1*(x8_8_5 - 1.0)*(x8_8_5 - 1.0) + 5.0d-1*(y8_8_5 - 1.0)*(y8_8_5 - 1.0) + 
	5.0d-1*(z8_8_5 - 1.0)*(z8_8_5 - 1.0) + 5.0d-1*(x9_8_5 - 1.0)*(x9_8_5 - 1.0) + 
	5.0d-1*(y9_8_5 - 1.0)*(y9_8_5 - 1.0) + 5.0d-1*(z9_8_5 - 1.0)*(z9_8_5 - 1.0) + 
	5.0d-1*(x1_9_5 - 1.0)*(x1_9_5 - 1.0) + 5.0d-1*(y1_9_5 - 1.0)*(y1_9_5 - 1.0) + 
	5.0d-1*(z1_9_5 - 1.0)*(z1_9_5 - 1.0) + 5.0d-1*(x2_9_5 - 1.0)*(x2_9_5 - 1.0) + 
	5.0d-1*(y2_9_5 - 1.0)*(y2_9_5 - 1.0) + 5.0d-1*(z2_9_5 - 1.0)*(z2_9_5 - 1.0) + 
	5.0d-1*(x3_9_5 - 1.0)*(x3_9_5 - 1.0) + 5.0d-1*(y3_9_5 - 1.0)*(y3_9_5 - 1.0) + 
	5.0d-1*(z3_9_5 - 1.0)*(z3_9_5 - 1.0) + 5.0d-1*(x4_9_5 - 1.0)*(x4_9_5 - 1.0) + 
	5.0d-1*(y4_9_5 - 1.0)*(y4_9_5 - 1.0) + 5.0d-1*(z4_9_5 - 1.0)*(z4_9_5 - 1.0) + 
	5.0d-1*(x5_9_5 - 1.0)*(x5_9_5 - 1.0) + 5.0d-1*(y5_9_5 - 1.0)*(y5_9_5 - 1.0) + 
	5.0d-1*(z5_9_5 - 1.0)*(z5_9_5 - 1.0) + 5.0d-1*(x6_9_5 - 1.0)*(x6_9_5 - 1.0) + 
	5.0d-1*(y6_9_5 - 1.0)*(y6_9_5 - 1.0) + 5.0d-1*(z6_9_5 - 1.0)*(z6_9_5 - 1.0) + 
	5.0d-1*(x7_9_5 - 1.0)*(x7_9_5 - 1.0) + 5.0d-1*(y7_9_5 - 1.0)*(y7_9_5 - 1.0) + 
	5.0d-1*(z7_9_5 - 1.0)*(z7_9_5 - 1.0) + 5.0d-1*(x8_9_5 - 1.0)*(x8_9_5 - 1.0) + 
	5.0d-1*(y8_9_5 - 1.0)*(y8_9_5 - 1.0) + 5.0d-1*(z8_9_5 - 1.0)*(z8_9_5 - 1.0) + 
	5.0d-1*(x9_9_5 - 1.0)*(x9_9_5 - 1.0) + 5.0d-1*(y9_9_5 - 1.0)*(y9_9_5 - 1.0) + 
	5.0d-1*(z9_9_5 - 1.0)*(z9_9_5 - 1.0) + 5.0d-1*(x1_1_6 - 1.0)*(x1_1_6 - 1.0) + 
	5.0d-1*(y1_1_6 - 1.0)*(y1_1_6 - 1.0) + 5.0d-1*(z1_1_6 - 1.0)*(z1_1_6 - 1.0) + 
	5.0d-1*(x2_1_6 - 1.0)*(x2_1_6 - 1.0) + 5.0d-1*(y2_1_6 - 1.0)*(y2_1_6 - 1.0) + 
	5.0d-1*(z2_1_6 - 1.0)*(z2_1_6 - 1.0) + 5.0d-1*(x3_1_6 - 1.0)*(x3_1_6 - 1.0) + 
	5.0d-1*(y3_1_6 - 1.0)*(y3_1_6 - 1.0) + 5.0d-1*(z3_1_6 - 1.0)*(z3_1_6 - 1.0) + 
	5.0d-1*(x4_1_6 - 1.0)*(x4_1_6 - 1.0) + 5.0d-1*(y4_1_6 - 1.0)*(y4_1_6 - 1.0) + 
	5.0d-1*(z4_1_6 - 1.0)*(z4_1_6 - 1.0) + 5.0d-1*(x5_1_6 - 1.0)*(x5_1_6 - 1.0) + 
	5.0d-1*(y5_1_6 - 1.0)*(y5_1_6 - 1.0) + 5.0d-1*(z5_1_6 - 1.0)*(z5_1_6 - 1.0) + 
	5.0d-1*(x6_1_6 - 1.0)*(x6_1_6 - 1.0) + 5.0d-1*(y6_1_6 - 1.0)*(y6_1_6 - 1.0) + 
	5.0d-1*(z6_1_6 - 1.0)*(z6_1_6 - 1.0) + 5.0d-1*(x7_1_6 - 1.0)*(x7_1_6 - 1.0) + 
	5.0d-1*(y7_1_6 - 1.0)*(y7_1_6 - 1.0) + 5.0d-1*(z7_1_6 - 1.0)*(z7_1_6 - 1.0) + 
	5.0d-1*(x8_1_6 - 1.0)*(x8_1_6 - 1.0) + 5.0d-1*(y8_1_6 - 1.0)*(y8_1_6 - 1.0) + 
	5.0d-1*(z8_1_6 - 1.0)*(z8_1_6 - 1.0) + 5.0d-1*(x9_1_6 - 1.0)*(x9_1_6 - 1.0) + 
	5.0d-1*(y9_1_6 - 1.0)*(y9_1_6 - 1.0) + 5.0d-1*(z9_1_6 - 1.0)*(z9_1_6 - 1.0) + 
	5.0d-1*(x1_2_6 - 1.0)*(x1_2_6 - 1.0) + 5.0d-1*(y1_2_6 - 1.0)*(y1_2_6 - 1.0) + 
	5.0d-1*(z1_2_6 - 1.0)*(z1_2_6 - 1.0) + 5.0d-1*(x2_2_6 - 1.0)*(x2_2_6 - 1.0) + 
	5.0d-1*(y2_2_6 - 1.0)*(y2_2_6 - 1.0) + 5.0d-1*(z2_2_6 - 1.0)*(z2_2_6 - 1.0) + 
	5.0d-1*(x3_2_6 - 1.0)*(x3_2_6 - 1.0) + 5.0d-1*(y3_2_6 - 1.0)*(y3_2_6 - 1.0) + 
	5.0d-1*(z3_2_6 - 1.0)*(z3_2_6 - 1.0) + 5.0d-1*(x4_2_6 - 1.0)*(x4_2_6 - 1.0) + 
	5.0d-1*(y4_2_6 - 1.0)*(y4_2_6 - 1.0) + 5.0d-1*(z4_2_6 - 1.0)*(z4_2_6 - 1.0) + 
	5.0d-1*(x5_2_6 - 1.0)*(x5_2_6 - 1.0) + 5.0d-1*(y5_2_6 - 1.0)*(y5_2_6 - 1.0) + 
	5.0d-1*(z5_2_6 - 1.0)*(z5_2_6 - 1.0) + 5.0d-1*(x6_2_6 - 1.0)*(x6_2_6 - 1.0) + 
	5.0d-1*(y6_2_6 - 1.0)*(y6_2_6 - 1.0) + 5.0d-1*(z6_2_6 - 1.0)*(z6_2_6 - 1.0) + 
	5.0d-1*(x7_2_6 - 1.0)*(x7_2_6 - 1.0) + 5.0d-1*(y7_2_6 - 1.0)*(y7_2_6 - 1.0) + 
	5.0d-1*(z7_2_6 - 1.0)*(z7_2_6 - 1.0) + 5.0d-1*(x8_2_6 - 1.0)*(x8_2_6 - 1.0) + 
	5.0d-1*(y8_2_6 - 1.0)*(y8_2_6 - 1.0) + 5.0d-1*(z8_2_6 - 1.0)*(z8_2_6 - 1.0) + 
	5.0d-1*(x9_2_6 - 1.0)*(x9_2_6 - 1.0) + 5.0d-1*(y9_2_6 - 1.0)*(y9_2_6 - 1.0) + 
	5.0d-1*(z9_2_6 - 1.0)*(z9_2_6 - 1.0) + 5.0d-1*(x1_3_6 - 1.0)*(x1_3_6 - 1.0) + 
	5.0d-1*(y1_3_6 - 1.0)*(y1_3_6 - 1.0) + 5.0d-1*(z1_3_6 - 1.0)*(z1_3_6 - 1.0) + 
	5.0d-1*(x2_3_6 - 1.0)*(x2_3_6 - 1.0) + 5.0d-1*(y2_3_6 - 1.0)*(y2_3_6 - 1.0) + 
	5.0d-1*(z2_3_6 - 1.0)*(z2_3_6 - 1.0) + 5.0d-1*(x3_3_6 - 1.0)*(x3_3_6 - 1.0) + 
	5.0d-1*(y3_3_6 - 1.0)*(y3_3_6 - 1.0) + 5.0d-1*(z3_3_6 - 1.0)*(z3_3_6 - 1.0) + 
	5.0d-1*(x4_3_6 - 1.0)*(x4_3_6 - 1.0) + 5.0d-1*(y4_3_6 - 1.0)*(y4_3_6 - 1.0) + 
	5.0d-1*(z4_3_6 - 1.0)*(z4_3_6 - 1.0) + 5.0d-1*(x5_3_6 - 1.0)*(x5_3_6 - 1.0) + 
	5.0d-1*(y5_3_6 - 1.0)*(y5_3_6 - 1.0) + 5.0d-1*(z5_3_6 - 1.0)*(z5_3_6 - 1.0) + 
	5.0d-1*(x6_3_6 - 1.0)*(x6_3_6 - 1.0) + 5.0d-1*(y6_3_6 - 1.0)*(y6_3_6 - 1.0) + 
	5.0d-1*(z6_3_6 - 1.0)*(z6_3_6 - 1.0) + 5.0d-1*(x7_3_6 - 1.0)*(x7_3_6 - 1.0) + 
	5.0d-1*(y7_3_6 - 1.0)*(y7_3_6 - 1.0) + 5.0d-1*(z7_3_6 - 1.0)*(z7_3_6 - 1.0) + 
	5.0d-1*(x8_3_6 - 1.0)*(x8_3_6 - 1.0) + 5.0d-1*(y8_3_6 - 1.0)*(y8_3_6 - 1.0) + 
	5.0d-1*(z8_3_6 - 1.0)*(z8_3_6 - 1.0) + 5.0d-1*(x9_3_6 - 1.0)*(x9_3_6 - 1.0) + 
	5.0d-1*(y9_3_6 - 1.0)*(y9_3_6 - 1.0) + 5.0d-1*(z9_3_6 - 1.0)*(z9_3_6 - 1.0) + 
	5.0d-1*(x1_4_6 - 1.0)*(x1_4_6 - 1.0) + 5.0d-1*(y1_4_6 - 1.0)*(y1_4_6 - 1.0) + 
	5.0d-1*(z1_4_6 - 1.0)*(z1_4_6 - 1.0) + 5.0d-1*(x2_4_6 - 1.0)*(x2_4_6 - 1.0) + 
	5.0d-1*(y2_4_6 - 1.0)*(y2_4_6 - 1.0) + 5.0d-1*(z2_4_6 - 1.0)*(z2_4_6 - 1.0) + 
	5.0d-1*(x3_4_6 - 1.0)*(x3_4_6 - 1.0) + 5.0d-1*(y3_4_6 - 1.0)*(y3_4_6 - 1.0) + 
	5.0d-1*(z3_4_6 - 1.0)*(z3_4_6 - 1.0) + 5.0d-1*(x4_4_6 - 1.0)*(x4_4_6 - 1.0) + 
	5.0d-1*(y4_4_6 - 1.0)*(y4_4_6 - 1.0) + 5.0d-1*(z4_4_6 - 1.0)*(z4_4_6 - 1.0) + 
	5.0d-1*(x5_4_6 - 1.0)*(x5_4_6 - 1.0) + 5.0d-1*(y5_4_6 - 1.0)*(y5_4_6 - 1.0) + 
	5.0d-1*(z5_4_6 - 1.0)*(z5_4_6 - 1.0) + 5.0d-1*(x6_4_6 - 1.0)*(x6_4_6 - 1.0) + 
	5.0d-1*(y6_4_6 - 1.0)*(y6_4_6 - 1.0) + 5.0d-1*(z6_4_6 - 1.0)*(z6_4_6 - 1.0) + 
	5.0d-1*(x7_4_6 - 1.0)*(x7_4_6 - 1.0) + 5.0d-1*(y7_4_6 - 1.0)*(y7_4_6 - 1.0) + 
	5.0d-1*(z7_4_6 - 1.0)*(z7_4_6 - 1.0) + 5.0d-1*(x8_4_6 - 1.0)*(x8_4_6 - 1.0) + 
	5.0d-1*(y8_4_6 - 1.0)*(y8_4_6 - 1.0) + 5.0d-1*(z8_4_6 - 1.0)*(z8_4_6 - 1.0) + 
	5.0d-1*(x9_4_6 - 1.0)*(x9_4_6 - 1.0) + 5.0d-1*(y9_4_6 - 1.0)*(y9_4_6 - 1.0) + 
	5.0d-1*(z9_4_6 - 1.0)*(z9_4_6 - 1.0) + 5.0d-1*(x1_5_6 - 1.0)*(x1_5_6 - 1.0) + 
	5.0d-1*(y1_5_6 - 1.0)*(y1_5_6 - 1.0) + 5.0d-1*(z1_5_6 - 1.0)*(z1_5_6 - 1.0) + 
	5.0d-1*(x2_5_6 - 1.0)*(x2_5_6 - 1.0) + 5.0d-1*(y2_5_6 - 1.0)*(y2_5_6 - 1.0) + 
	5.0d-1*(z2_5_6 - 1.0)*(z2_5_6 - 1.0) + 5.0d-1*(x3_5_6 - 1.0)*(x3_5_6 - 1.0) + 
	5.0d-1*(y3_5_6 - 1.0)*(y3_5_6 - 1.0) + 5.0d-1*(z3_5_6 - 1.0)*(z3_5_6 - 1.0) + 
	5.0d-1*(x4_5_6 - 1.0)*(x4_5_6 - 1.0) + 5.0d-1*(y4_5_6 - 1.0)*(y4_5_6 - 1.0) + 
	5.0d-1*(z4_5_6 - 1.0)*(z4_5_6 - 1.0) + 5.0d-1*(x5_5_6 - 1.0)*(x5_5_6 - 1.0) + 
	5.0d-1*(y5_5_6 - 1.0)*(y5_5_6 - 1.0) + 5.0d-1*(z5_5_6 - 1.0)*(z5_5_6 - 1.0) + 
	5.0d-1*(x6_5_6 - 1.0)*(x6_5_6 - 1.0) + 5.0d-1*(y6_5_6 - 1.0)*(y6_5_6 - 1.0) + 
	5.0d-1*(z6_5_6 - 1.0)*(z6_5_6 - 1.0) + 5.0d-1*(x7_5_6 - 1.0)*(x7_5_6 - 1.0) + 
	5.0d-1*(y7_5_6 - 1.0)*(y7_5_6 - 1.0) + 5.0d-1*(z7_5_6 - 1.0)*(z7_5_6 - 1.0) + 
	5.0d-1*(x8_5_6 - 1.0)*(x8_5_6 - 1.0) + 5.0d-1*(y8_5_6 - 1.0)*(y8_5_6 - 1.0) + 
	5.0d-1*(z8_5_6 - 1.0)*(z8_5_6 - 1.0) + 5.0d-1*(x9_5_6 - 1.0)*(x9_5_6 - 1.0) + 
	5.0d-1*(y9_5_6 - 1.0)*(y9_5_6 - 1.0) + 5.0d-1*(z9_5_6 - 1.0)*(z9_5_6 - 1.0) + 
	5.0d-1*(x1_6_6 - 1.0)*(x1_6_6 - 1.0) + 5.0d-1*(y1_6_6 - 1.0)*(y1_6_6 - 1.0) + 
	5.0d-1*(z1_6_6 - 1.0)*(z1_6_6 - 1.0) + 5.0d-1*(x2_6_6 - 1.0)*(x2_6_6 - 1.0) + 
	5.0d-1*(y2_6_6 - 1.0)*(y2_6_6 - 1.0) + 5.0d-1*(z2_6_6 - 1.0)*(z2_6_6 - 1.0) + 
	5.0d-1*(x3_6_6 - 1.0)*(x3_6_6 - 1.0) + 5.0d-1*(y3_6_6 - 1.0)*(y3_6_6 - 1.0) + 
	5.0d-1*(z3_6_6 - 1.0)*(z3_6_6 - 1.0) + 5.0d-1*(x4_6_6 - 1.0)*(x4_6_6 - 1.0) + 
	5.0d-1*(y4_6_6 - 1.0)*(y4_6_6 - 1.0) + 5.0d-1*(z4_6_6 - 1.0)*(z4_6_6 - 1.0) + 
	5.0d-1*(x5_6_6 - 1.0)*(x5_6_6 - 1.0) + 5.0d-1*(y5_6_6 - 1.0)*(y5_6_6 - 1.0) + 
	5.0d-1*(z5_6_6 - 1.0)*(z5_6_6 - 1.0) + 5.0d-1*(x6_6_6 - 1.0)*(x6_6_6 - 1.0) + 
	5.0d-1*(y6_6_6 - 1.0)*(y6_6_6 - 1.0) + 5.0d-1*(z6_6_6 - 1.0)*(z6_6_6 - 1.0) + 
	5.0d-1*(x7_6_6 - 1.0)*(x7_6_6 - 1.0) + 5.0d-1*(y7_6_6 - 1.0)*(y7_6_6 - 1.0) + 
	5.0d-1*(z7_6_6 - 1.0)*(z7_6_6 - 1.0) + 5.0d-1*(x8_6_6 - 1.0)*(x8_6_6 - 1.0) + 
	5.0d-1*(y8_6_6 - 1.0)*(y8_6_6 - 1.0) + 5.0d-1*(z8_6_6 - 1.0)*(z8_6_6 - 1.0) + 
	5.0d-1*(x9_6_6 - 1.0)*(x9_6_6 - 1.0) + 5.0d-1*(y9_6_6 - 1.0)*(y9_6_6 - 1.0) + 
	5.0d-1*(z9_6_6 - 1.0)*(z9_6_6 - 1.0) + 5.0d-1*(x1_7_6 - 1.0)*(x1_7_6 - 1.0) + 
	5.0d-1*(y1_7_6 - 1.0)*(y1_7_6 - 1.0) + 5.0d-1*(z1_7_6 - 1.0)*(z1_7_6 - 1.0) + 
	5.0d-1*(x2_7_6 - 1.0)*(x2_7_6 - 1.0) + 5.0d-1*(y2_7_6 - 1.0)*(y2_7_6 - 1.0) + 
	5.0d-1*(z2_7_6 - 1.0)*(z2_7_6 - 1.0) + 5.0d-1*(x3_7_6 - 1.0)*(x3_7_6 - 1.0) + 
	5.0d-1*(y3_7_6 - 1.0)*(y3_7_6 - 1.0) + 5.0d-1*(z3_7_6 - 1.0)*(z3_7_6 - 1.0) + 
	5.0d-1*(x4_7_6 - 1.0)*(x4_7_6 - 1.0) + 5.0d-1*(y4_7_6 - 1.0)*(y4_7_6 - 1.0) + 
	5.0d-1*(z4_7_6 - 1.0)*(z4_7_6 - 1.0) + 5.0d-1*(x5_7_6 - 1.0)*(x5_7_6 - 1.0) + 
	5.0d-1*(y5_7_6 - 1.0)*(y5_7_6 - 1.0) + 5.0d-1*(z5_7_6 - 1.0)*(z5_7_6 - 1.0) + 
	5.0d-1*(x6_7_6 - 1.0)*(x6_7_6 - 1.0) + 5.0d-1*(y6_7_6 - 1.0)*(y6_7_6 - 1.0) + 
	5.0d-1*(z6_7_6 - 1.0)*(z6_7_6 - 1.0) + 5.0d-1*(x7_7_6 - 1.0)*(x7_7_6 - 1.0) + 
	5.0d-1*(y7_7_6 - 1.0)*(y7_7_6 - 1.0) + 5.0d-1*(z7_7_6 - 1.0)*(z7_7_6 - 1.0) + 
	5.0d-1*(x8_7_6 - 1.0)*(x8_7_6 - 1.0) + 5.0d-1*(y8_7_6 - 1.0)*(y8_7_6 - 1.0) + 
	5.0d-1*(z8_7_6 - 1.0)*(z8_7_6 - 1.0) + 5.0d-1*(x9_7_6 - 1.0)*(x9_7_6 - 1.0) + 
	5.0d-1*(y9_7_6 - 1.0)*(y9_7_6 - 1.0) + 5.0d-1*(z9_7_6 - 1.0)*(z9_7_6 - 1.0) + 
	5.0d-1*(x1_8_6 - 1.0)*(x1_8_6 - 1.0) + 5.0d-1*(y1_8_6 - 1.0)*(y1_8_6 - 1.0) + 
	5.0d-1*(z1_8_6 - 1.0)*(z1_8_6 - 1.0) + 5.0d-1*(x2_8_6 - 1.0)*(x2_8_6 - 1.0) + 
	5.0d-1*(y2_8_6 - 1.0)*(y2_8_6 - 1.0) + 5.0d-1*(z2_8_6 - 1.0)*(z2_8_6 - 1.0) + 
	5.0d-1*(x3_8_6 - 1.0)*(x3_8_6 - 1.0) + 5.0d-1*(y3_8_6 - 1.0)*(y3_8_6 - 1.0) + 
	5.0d-1*(z3_8_6 - 1.0)*(z3_8_6 - 1.0) + 5.0d-1*(x4_8_6 - 1.0)*(x4_8_6 - 1.0) + 
	5.0d-1*(y4_8_6 - 1.0)*(y4_8_6 - 1.0) + 5.0d-1*(z4_8_6 - 1.0)*(z4_8_6 - 1.0) + 
	5.0d-1*(x5_8_6 - 1.0)*(x5_8_6 - 1.0) + 5.0d-1*(y5_8_6 - 1.0)*(y5_8_6 - 1.0) + 
	5.0d-1*(z5_8_6 - 1.0)*(z5_8_6 - 1.0) + 5.0d-1*(x6_8_6 - 1.0)*(x6_8_6 - 1.0) + 
	5.0d-1*(y6_8_6 - 1.0)*(y6_8_6 - 1.0) + 5.0d-1*(z6_8_6 - 1.0)*(z6_8_6 - 1.0) + 
	5.0d-1*(x7_8_6 - 1.0)*(x7_8_6 - 1.0) + 5.0d-1*(y7_8_6 - 1.0)*(y7_8_6 - 1.0) + 
	5.0d-1*(z7_8_6 - 1.0)*(z7_8_6 - 1.0) + 5.0d-1*(x8_8_6 - 1.0)*(x8_8_6 - 1.0) + 
	5.0d-1*(y8_8_6 - 1.0)*(y8_8_6 - 1.0) + 5.0d-1*(z8_8_6 - 1.0)*(z8_8_6 - 1.0) + 
	5.0d-1*(x9_8_6 - 1.0)*(x9_8_6 - 1.0) + 5.0d-1*(y9_8_6 - 1.0)*(y9_8_6 - 1.0) + 
	5.0d-1*(z9_8_6 - 1.0)*(z9_8_6 - 1.0) + 5.0d-1*(x1_9_6 - 1.0)*(x1_9_6 - 1.0) + 
	5.0d-1*(y1_9_6 - 1.0)*(y1_9_6 - 1.0) + 5.0d-1*(z1_9_6 - 1.0)*(z1_9_6 - 1.0) + 
	5.0d-1*(x2_9_6 - 1.0)*(x2_9_6 - 1.0) + 5.0d-1*(y2_9_6 - 1.0)*(y2_9_6 - 1.0) + 
	5.0d-1*(z2_9_6 - 1.0)*(z2_9_6 - 1.0) + 5.0d-1*(x3_9_6 - 1.0)*(x3_9_6 - 1.0) + 
	5.0d-1*(y3_9_6 - 1.0)*(y3_9_6 - 1.0) + 5.0d-1*(z3_9_6 - 1.0)*(z3_9_6 - 1.0) + 
	5.0d-1*(x4_9_6 - 1.0)*(x4_9_6 - 1.0) + 5.0d-1*(y4_9_6 - 1.0)*(y4_9_6 - 1.0) + 
	5.0d-1*(z4_9_6 - 1.0)*(z4_9_6 - 1.0) + 5.0d-1*(x5_9_6 - 1.0)*(x5_9_6 - 1.0) + 
	5.0d-1*(y5_9_6 - 1.0)*(y5_9_6 - 1.0) + 5.0d-1*(z5_9_6 - 1.0)*(z5_9_6 - 1.0) + 
	5.0d-1*(x6_9_6 - 1.0)*(x6_9_6 - 1.0) + 5.0d-1*(y6_9_6 - 1.0)*(y6_9_6 - 1.0) + 
	5.0d-1*(z6_9_6 - 1.0)*(z6_9_6 - 1.0) + 5.0d-1*(x7_9_6 - 1.0)*(x7_9_6 - 1.0) + 
	5.0d-1*(y7_9_6 - 1.0)*(y7_9_6 - 1.0) + 5.0d-1*(z7_9_6 - 1.0)*(z7_9_6 - 1.0) + 
	5.0d-1*(x8_9_6 - 1.0)*(x8_9_6 - 1.0) + 5.0d-1*(y8_9_6 - 1.0)*(y8_9_6 - 1.0) + 
	5.0d-1*(z8_9_6 - 1.0)*(z8_9_6 - 1.0) + 5.0d-1*(x9_9_6 - 1.0)*(x9_9_6 - 1.0) + 
	5.0d-1*(y9_9_6 - 1.0)*(y9_9_6 - 1.0) + 5.0d-1*(z9_9_6 - 1.0)*(z9_9_6 - 1.0) + 
	5.0d-1*(x1_1_7 - 1.0)*(x1_1_7 - 1.0) + 5.0d-1*(y1_1_7 - 1.0)*(y1_1_7 - 1.0) + 
	5.0d-1*(z1_1_7 - 1.0)*(z1_1_7 - 1.0) + 5.0d-1*(x2_1_7 - 1.0)*(x2_1_7 - 1.0) + 
	5.0d-1*(y2_1_7 - 1.0)*(y2_1_7 - 1.0) + 5.0d-1*(z2_1_7 - 1.0)*(z2_1_7 - 1.0) + 
	5.0d-1*(x3_1_7 - 1.0)*(x3_1_7 - 1.0) + 5.0d-1*(y3_1_7 - 1.0)*(y3_1_7 - 1.0) + 
	5.0d-1*(z3_1_7 - 1.0)*(z3_1_7 - 1.0) + 5.0d-1*(x4_1_7 - 1.0)*(x4_1_7 - 1.0) + 
	5.0d-1*(y4_1_7 - 1.0)*(y4_1_7 - 1.0) + 5.0d-1*(z4_1_7 - 1.0)*(z4_1_7 - 1.0) + 
	5.0d-1*(x5_1_7 - 1.0)*(x5_1_7 - 1.0) + 5.0d-1*(y5_1_7 - 1.0)*(y5_1_7 - 1.0) + 
	5.0d-1*(z5_1_7 - 1.0)*(z5_1_7 - 1.0) + 5.0d-1*(x6_1_7 - 1.0)*(x6_1_7 - 1.0) + 
	5.0d-1*(y6_1_7 - 1.0)*(y6_1_7 - 1.0) + 5.0d-1*(z6_1_7 - 1.0)*(z6_1_7 - 1.0) + 
	5.0d-1*(x7_1_7 - 1.0)*(x7_1_7 - 1.0) + 5.0d-1*(y7_1_7 - 1.0)*(y7_1_7 - 1.0) + 
	5.0d-1*(z7_1_7 - 1.0)*(z7_1_7 - 1.0) + 5.0d-1*(x8_1_7 - 1.0)*(x8_1_7 - 1.0) + 
	5.0d-1*(y8_1_7 - 1.0)*(y8_1_7 - 1.0) + 5.0d-1*(z8_1_7 - 1.0)*(z8_1_7 - 1.0) + 
	5.0d-1*(x9_1_7 - 1.0)*(x9_1_7 - 1.0) + 5.0d-1*(y9_1_7 - 1.0)*(y9_1_7 - 1.0) + 
	5.0d-1*(z9_1_7 - 1.0)*(z9_1_7 - 1.0) + 5.0d-1*(x1_2_7 - 1.0)*(x1_2_7 - 1.0) + 
	5.0d-1*(y1_2_7 - 1.0)*(y1_2_7 - 1.0) + 5.0d-1*(z1_2_7 - 1.0)*(z1_2_7 - 1.0) + 
	5.0d-1*(x2_2_7 - 1.0)*(x2_2_7 - 1.0) + 5.0d-1*(y2_2_7 - 1.0)*(y2_2_7 - 1.0) + 
	5.0d-1*(z2_2_7 - 1.0)*(z2_2_7 - 1.0) + 5.0d-1*(x3_2_7 - 1.0)*(x3_2_7 - 1.0) + 
	5.0d-1*(y3_2_7 - 1.0)*(y3_2_7 - 1.0) + 5.0d-1*(z3_2_7 - 1.0)*(z3_2_7 - 1.0) + 
	5.0d-1*(x4_2_7 - 1.0)*(x4_2_7 - 1.0) + 5.0d-1*(y4_2_7 - 1.0)*(y4_2_7 - 1.0) + 
	5.0d-1*(z4_2_7 - 1.0)*(z4_2_7 - 1.0) + 5.0d-1*(x5_2_7 - 1.0)*(x5_2_7 - 1.0) + 
	5.0d-1*(y5_2_7 - 1.0)*(y5_2_7 - 1.0) + 5.0d-1*(z5_2_7 - 1.0)*(z5_2_7 - 1.0) + 
	5.0d-1*(x6_2_7 - 1.0)*(x6_2_7 - 1.0) + 5.0d-1*(y6_2_7 - 1.0)*(y6_2_7 - 1.0) + 
	5.0d-1*(z6_2_7 - 1.0)*(z6_2_7 - 1.0) + 5.0d-1*(x7_2_7 - 1.0)*(x7_2_7 - 1.0) + 
	5.0d-1*(y7_2_7 - 1.0)*(y7_2_7 - 1.0) + 5.0d-1*(z7_2_7 - 1.0)*(z7_2_7 - 1.0) + 
	5.0d-1*(x8_2_7 - 1.0)*(x8_2_7 - 1.0) + 5.0d-1*(y8_2_7 - 1.0)*(y8_2_7 - 1.0) + 
	5.0d-1*(z8_2_7 - 1.0)*(z8_2_7 - 1.0) + 5.0d-1*(x9_2_7 - 1.0)*(x9_2_7 - 1.0) + 
	5.0d-1*(y9_2_7 - 1.0)*(y9_2_7 - 1.0) + 5.0d-1*(z9_2_7 - 1.0)*(z9_2_7 - 1.0) + 
	5.0d-1*(x1_3_7 - 1.0)*(x1_3_7 - 1.0) + 5.0d-1*(y1_3_7 - 1.0)*(y1_3_7 - 1.0) + 
	5.0d-1*(z1_3_7 - 1.0)*(z1_3_7 - 1.0) + 5.0d-1*(x2_3_7 - 1.0)*(x2_3_7 - 1.0) + 
	5.0d-1*(y2_3_7 - 1.0)*(y2_3_7 - 1.0) + 5.0d-1*(z2_3_7 - 1.0)*(z2_3_7 - 1.0) + 
	5.0d-1*(x3_3_7 - 1.0)*(x3_3_7 - 1.0) + 5.0d-1*(y3_3_7 - 1.0)*(y3_3_7 - 1.0) + 
	5.0d-1*(z3_3_7 - 1.0)*(z3_3_7 - 1.0) + 5.0d-1*(x4_3_7 - 1.0)*(x4_3_7 - 1.0) + 
	5.0d-1*(y4_3_7 - 1.0)*(y4_3_7 - 1.0) + 5.0d-1*(z4_3_7 - 1.0)*(z4_3_7 - 1.0) + 
	5.0d-1*(x5_3_7 - 1.0)*(x5_3_7 - 1.0) + 5.0d-1*(y5_3_7 - 1.0)*(y5_3_7 - 1.0) + 
	5.0d-1*(z5_3_7 - 1.0)*(z5_3_7 - 1.0) + 5.0d-1*(x6_3_7 - 1.0)*(x6_3_7 - 1.0) + 
	5.0d-1*(y6_3_7 - 1.0)*(y6_3_7 - 1.0) + 5.0d-1*(z6_3_7 - 1.0)*(z6_3_7 - 1.0) + 
	5.0d-1*(x7_3_7 - 1.0)*(x7_3_7 - 1.0) + 5.0d-1*(y7_3_7 - 1.0)*(y7_3_7 - 1.0) + 
	5.0d-1*(z7_3_7 - 1.0)*(z7_3_7 - 1.0) + 5.0d-1*(x8_3_7 - 1.0)*(x8_3_7 - 1.0) + 
	5.0d-1*(y8_3_7 - 1.0)*(y8_3_7 - 1.0) + 5.0d-1*(z8_3_7 - 1.0)*(z8_3_7 - 1.0) + 
	5.0d-1*(x9_3_7 - 1.0)*(x9_3_7 - 1.0) + 5.0d-1*(y9_3_7 - 1.0)*(y9_3_7 - 1.0) + 
	5.0d-1*(z9_3_7 - 1.0)*(z9_3_7 - 1.0) + 5.0d-1*(x1_4_7 - 1.0)*(x1_4_7 - 1.0) + 
	5.0d-1*(y1_4_7 - 1.0)*(y1_4_7 - 1.0) + 5.0d-1*(z1_4_7 - 1.0)*(z1_4_7 - 1.0) + 
	5.0d-1*(x2_4_7 - 1.0)*(x2_4_7 - 1.0) + 5.0d-1*(y2_4_7 - 1.0)*(y2_4_7 - 1.0) + 
	5.0d-1*(z2_4_7 - 1.0)*(z2_4_7 - 1.0) + 5.0d-1*(x3_4_7 - 1.0)*(x3_4_7 - 1.0) + 
	5.0d-1*(y3_4_7 - 1.0)*(y3_4_7 - 1.0) + 5.0d-1*(z3_4_7 - 1.0)*(z3_4_7 - 1.0) + 
	5.0d-1*(x4_4_7 - 1.0)*(x4_4_7 - 1.0) + 5.0d-1*(y4_4_7 - 1.0)*(y4_4_7 - 1.0) + 
	5.0d-1*(z4_4_7 - 1.0)*(z4_4_7 - 1.0) + 5.0d-1*(x5_4_7 - 1.0)*(x5_4_7 - 1.0) + 
	5.0d-1*(y5_4_7 - 1.0)*(y5_4_7 - 1.0) + 5.0d-1*(z5_4_7 - 1.0)*(z5_4_7 - 1.0) + 
	5.0d-1*(x6_4_7 - 1.0)*(x6_4_7 - 1.0) + 5.0d-1*(y6_4_7 - 1.0)*(y6_4_7 - 1.0) + 
	5.0d-1*(z6_4_7 - 1.0)*(z6_4_7 - 1.0) + 5.0d-1*(x7_4_7 - 1.0)*(x7_4_7 - 1.0) + 
	5.0d-1*(y7_4_7 - 1.0)*(y7_4_7 - 1.0) + 5.0d-1*(z7_4_7 - 1.0)*(z7_4_7 - 1.0) + 
	5.0d-1*(x8_4_7 - 1.0)*(x8_4_7 - 1.0) + 5.0d-1*(y8_4_7 - 1.0)*(y8_4_7 - 1.0) + 
	5.0d-1*(z8_4_7 - 1.0)*(z8_4_7 - 1.0) + 5.0d-1*(x9_4_7 - 1.0)*(x9_4_7 - 1.0) + 
	5.0d-1*(y9_4_7 - 1.0)*(y9_4_7 - 1.0) + 5.0d-1*(z9_4_7 - 1.0)*(z9_4_7 - 1.0) + 
	5.0d-1*(x1_5_7 - 1.0)*(x1_5_7 - 1.0) + 5.0d-1*(y1_5_7 - 1.0)*(y1_5_7 - 1.0) + 
	5.0d-1*(z1_5_7 - 1.0)*(z1_5_7 - 1.0) + 5.0d-1*(x2_5_7 - 1.0)*(x2_5_7 - 1.0) + 
	5.0d-1*(y2_5_7 - 1.0)*(y2_5_7 - 1.0) + 5.0d-1*(z2_5_7 - 1.0)*(z2_5_7 - 1.0) + 
	5.0d-1*(x3_5_7 - 1.0)*(x3_5_7 - 1.0) + 5.0d-1*(y3_5_7 - 1.0)*(y3_5_7 - 1.0) + 
	5.0d-1*(z3_5_7 - 1.0)*(z3_5_7 - 1.0) + 5.0d-1*(x4_5_7 - 1.0)*(x4_5_7 - 1.0) + 
	5.0d-1*(y4_5_7 - 1.0)*(y4_5_7 - 1.0) + 5.0d-1*(z4_5_7 - 1.0)*(z4_5_7 - 1.0) + 
	5.0d-1*(x5_5_7 - 1.0)*(x5_5_7 - 1.0) + 5.0d-1*(y5_5_7 - 1.0)*(y5_5_7 - 1.0) + 
	5.0d-1*(z5_5_7 - 1.0)*(z5_5_7 - 1.0) + 5.0d-1*(x6_5_7 - 1.0)*(x6_5_7 - 1.0) + 
	5.0d-1*(y6_5_7 - 1.0)*(y6_5_7 - 1.0) + 5.0d-1*(z6_5_7 - 1.0)*(z6_5_7 - 1.0) + 
	5.0d-1*(x7_5_7 - 1.0)*(x7_5_7 - 1.0) + 5.0d-1*(y7_5_7 - 1.0)*(y7_5_7 - 1.0) + 
	5.0d-1*(z7_5_7 - 1.0)*(z7_5_7 - 1.0) + 5.0d-1*(x8_5_7 - 1.0)*(x8_5_7 - 1.0) + 
	5.0d-1*(y8_5_7 - 1.0)*(y8_5_7 - 1.0) + 5.0d-1*(z8_5_7 - 1.0)*(z8_5_7 - 1.0) + 
	5.0d-1*(x9_5_7 - 1.0)*(x9_5_7 - 1.0) + 5.0d-1*(y9_5_7 - 1.0)*(y9_5_7 - 1.0) + 
	5.0d-1*(z9_5_7 - 1.0)*(z9_5_7 - 1.0) + 5.0d-1*(x1_6_7 - 1.0)*(x1_6_7 - 1.0) + 
	5.0d-1*(y1_6_7 - 1.0)*(y1_6_7 - 1.0) + 5.0d-1*(z1_6_7 - 1.0)*(z1_6_7 - 1.0) + 
	5.0d-1*(x2_6_7 - 1.0)*(x2_6_7 - 1.0) + 5.0d-1*(y2_6_7 - 1.0)*(y2_6_7 - 1.0) + 
	5.0d-1*(z2_6_7 - 1.0)*(z2_6_7 - 1.0) + 5.0d-1*(x3_6_7 - 1.0)*(x3_6_7 - 1.0) + 
	5.0d-1*(y3_6_7 - 1.0)*(y3_6_7 - 1.0) + 5.0d-1*(z3_6_7 - 1.0)*(z3_6_7 - 1.0) + 
	5.0d-1*(x4_6_7 - 1.0)*(x4_6_7 - 1.0) + 5.0d-1*(y4_6_7 - 1.0)*(y4_6_7 - 1.0) + 
	5.0d-1*(z4_6_7 - 1.0)*(z4_6_7 - 1.0) + 5.0d-1*(x5_6_7 - 1.0)*(x5_6_7 - 1.0) + 
	5.0d-1*(y5_6_7 - 1.0)*(y5_6_7 - 1.0) + 5.0d-1*(z5_6_7 - 1.0)*(z5_6_7 - 1.0) + 
	5.0d-1*(x6_6_7 - 1.0)*(x6_6_7 - 1.0) + 5.0d-1*(y6_6_7 - 1.0)*(y6_6_7 - 1.0) + 
	5.0d-1*(z6_6_7 - 1.0)*(z6_6_7 - 1.0) + 5.0d-1*(x7_6_7 - 1.0)*(x7_6_7 - 1.0) + 
	5.0d-1*(y7_6_7 - 1.0)*(y7_6_7 - 1.0) + 5.0d-1*(z7_6_7 - 1.0)*(z7_6_7 - 1.0) + 
	5.0d-1*(x8_6_7 - 1.0)*(x8_6_7 - 1.0) + 5.0d-1*(y8_6_7 - 1.0)*(y8_6_7 - 1.0) + 
	5.0d-1*(z8_6_7 - 1.0)*(z8_6_7 - 1.0) + 5.0d-1*(x9_6_7 - 1.0)*(x9_6_7 - 1.0) + 
	5.0d-1*(y9_6_7 - 1.0)*(y9_6_7 - 1.0) + 5.0d-1*(z9_6_7 - 1.0)*(z9_6_7 - 1.0) + 
	5.0d-1*(x1_7_7 - 1.0)*(x1_7_7 - 1.0) + 5.0d-1*(y1_7_7 - 1.0)*(y1_7_7 - 1.0) + 
	5.0d-1*(z1_7_7 - 1.0)*(z1_7_7 - 1.0) + 5.0d-1*(x2_7_7 - 1.0)*(x2_7_7 - 1.0) + 
	5.0d-1*(y2_7_7 - 1.0)*(y2_7_7 - 1.0) + 5.0d-1*(z2_7_7 - 1.0)*(z2_7_7 - 1.0) + 
	5.0d-1*(x3_7_7 - 1.0)*(x3_7_7 - 1.0) + 5.0d-1*(y3_7_7 - 1.0)*(y3_7_7 - 1.0) + 
	5.0d-1*(z3_7_7 - 1.0)*(z3_7_7 - 1.0) + 5.0d-1*(x4_7_7 - 1.0)*(x4_7_7 - 1.0) + 
	5.0d-1*(y4_7_7 - 1.0)*(y4_7_7 - 1.0) + 5.0d-1*(z4_7_7 - 1.0)*(z4_7_7 - 1.0) + 
	5.0d-1*(x5_7_7 - 1.0)*(x5_7_7 - 1.0) + 5.0d-1*(y5_7_7 - 1.0)*(y5_7_7 - 1.0) + 
	5.0d-1*(z5_7_7 - 1.0)*(z5_7_7 - 1.0) + 5.0d-1*(x6_7_7 - 1.0)*(x6_7_7 - 1.0) + 
	5.0d-1*(y6_7_7 - 1.0)*(y6_7_7 - 1.0) + 5.0d-1*(z6_7_7 - 1.0)*(z6_7_7 - 1.0) + 
	5.0d-1*(x7_7_7 - 1.0)*(x7_7_7 - 1.0) + 5.0d-1*(y7_7_7 - 1.0)*(y7_7_7 - 1.0) + 
	5.0d-1*(z7_7_7 - 1.0)*(z7_7_7 - 1.0) + 5.0d-1*(x8_7_7 - 1.0)*(x8_7_7 - 1.0) + 
	5.0d-1*(y8_7_7 - 1.0)*(y8_7_7 - 1.0) + 5.0d-1*(z8_7_7 - 1.0)*(z8_7_7 - 1.0) + 
	5.0d-1*(x9_7_7 - 1.0)*(x9_7_7 - 1.0) + 5.0d-1*(y9_7_7 - 1.0)*(y9_7_7 - 1.0) + 
	5.0d-1*(z9_7_7 - 1.0)*(z9_7_7 - 1.0) + 5.0d-1*(x1_8_7 - 1.0)*(x1_8_7 - 1.0) + 
	5.0d-1*(y1_8_7 - 1.0)*(y1_8_7 - 1.0) + 5.0d-1*(z1_8_7 - 1.0)*(z1_8_7 - 1.0) + 
	5.0d-1*(x2_8_7 - 1.0)*(x2_8_7 - 1.0) + 5.0d-1*(y2_8_7 - 1.0)*(y2_8_7 - 1.0) + 
	5.0d-1*(z2_8_7 - 1.0)*(z2_8_7 - 1.0) + 5.0d-1*(x3_8_7 - 1.0)*(x3_8_7 - 1.0) + 
	5.0d-1*(y3_8_7 - 1.0)*(y3_8_7 - 1.0) + 5.0d-1*(z3_8_7 - 1.0)*(z3_8_7 - 1.0) + 
	5.0d-1*(x4_8_7 - 1.0)*(x4_8_7 - 1.0) + 5.0d-1*(y4_8_7 - 1.0)*(y4_8_7 - 1.0) + 
	5.0d-1*(z4_8_7 - 1.0)*(z4_8_7 - 1.0) + 5.0d-1*(x5_8_7 - 1.0)*(x5_8_7 - 1.0) + 
	5.0d-1*(y5_8_7 - 1.0)*(y5_8_7 - 1.0) + 5.0d-1*(z5_8_7 - 1.0)*(z5_8_7 - 1.0) + 
	5.0d-1*(x6_8_7 - 1.0)*(x6_8_7 - 1.0) + 5.0d-1*(y6_8_7 - 1.0)*(y6_8_7 - 1.0) + 
	5.0d-1*(z6_8_7 - 1.0)*(z6_8_7 - 1.0) + 5.0d-1*(x7_8_7 - 1.0)*(x7_8_7 - 1.0) + 
	5.0d-1*(y7_8_7 - 1.0)*(y7_8_7 - 1.0) + 5.0d-1*(z7_8_7 - 1.0)*(z7_8_7 - 1.0) + 
	5.0d-1*(x8_8_7 - 1.0)*(x8_8_7 - 1.0) + 5.0d-1*(y8_8_7 - 1.0)*(y8_8_7 - 1.0) + 
	5.0d-1*(z8_8_7 - 1.0)*(z8_8_7 - 1.0) + 5.0d-1*(x9_8_7 - 1.0)*(x9_8_7 - 1.0) + 
	5.0d-1*(y9_8_7 - 1.0)*(y9_8_7 - 1.0) + 5.0d-1*(z9_8_7 - 1.0)*(z9_8_7 - 1.0) + 
	5.0d-1*(x1_9_7 - 1.0)*(x1_9_7 - 1.0) + 5.0d-1*(y1_9_7 - 1.0)*(y1_9_7 - 1.0) + 
	5.0d-1*(z1_9_7 - 1.0)*(z1_9_7 - 1.0) + 5.0d-1*(x2_9_7 - 1.0)*(x2_9_7 - 1.0) + 
	5.0d-1*(y2_9_7 - 1.0)*(y2_9_7 - 1.0) + 5.0d-1*(z2_9_7 - 1.0)*(z2_9_7 - 1.0) + 
	5.0d-1*(x3_9_7 - 1.0)*(x3_9_7 - 1.0) + 5.0d-1*(y3_9_7 - 1.0)*(y3_9_7 - 1.0) + 
	5.0d-1*(z3_9_7 - 1.0)*(z3_9_7 - 1.0) + 5.0d-1*(x4_9_7 - 1.0)*(x4_9_7 - 1.0) + 
	5.0d-1*(y4_9_7 - 1.0)*(y4_9_7 - 1.0) + 5.0d-1*(z4_9_7 - 1.0)*(z4_9_7 - 1.0) + 
	5.0d-1*(x5_9_7 - 1.0)*(x5_9_7 - 1.0) + 5.0d-1*(y5_9_7 - 1.0)*(y5_9_7 - 1.0) + 
	5.0d-1*(z5_9_7 - 1.0)*(z5_9_7 - 1.0) + 5.0d-1*(x6_9_7 - 1.0)*(x6_9_7 - 1.0) + 
	5.0d-1*(y6_9_7 - 1.0)*(y6_9_7 - 1.0) + 5.0d-1*(z6_9_7 - 1.0)*(z6_9_7 - 1.0) + 
	5.0d-1*(x7_9_7 - 1.0)*(x7_9_7 - 1.0) + 5.0d-1*(y7_9_7 - 1.0)*(y7_9_7 - 1.0) + 
	5.0d-1*(z7_9_7 - 1.0)*(z7_9_7 - 1.0) + 5.0d-1*(x8_9_7 - 1.0)*(x8_9_7 - 1.0) + 
	5.0d-1*(y8_9_7 - 1.0)*(y8_9_7 - 1.0) + 5.0d-1*(z8_9_7 - 1.0)*(z8_9_7 - 1.0) + 
	5.0d-1*(x9_9_7 - 1.0)*(x9_9_7 - 1.0) + 5.0d-1*(y9_9_7 - 1.0)*(y9_9_7 - 1.0) + 
	5.0d-1*(z9_9_7 - 1.0)*(z9_9_7 - 1.0) + 5.0d-1*(x1_1_8 - 1.0)*(x1_1_8 - 1.0) + 
	5.0d-1*(y1_1_8 - 1.0)*(y1_1_8 - 1.0) + 5.0d-1*(z1_1_8 - 1.0)*(z1_1_8 - 1.0) + 
	5.0d-1*(x2_1_8 - 1.0)*(x2_1_8 - 1.0) + 5.0d-1*(y2_1_8 - 1.0)*(y2_1_8 - 1.0) + 
	5.0d-1*(z2_1_8 - 1.0)*(z2_1_8 - 1.0) + 5.0d-1*(x3_1_8 - 1.0)*(x3_1_8 - 1.0) + 
	5.0d-1*(y3_1_8 - 1.0)*(y3_1_8 - 1.0) + 5.0d-1*(z3_1_8 - 1.0)*(z3_1_8 - 1.0) + 
	5.0d-1*(x4_1_8 - 1.0)*(x4_1_8 - 1.0) + 5.0d-1*(y4_1_8 - 1.0)*(y4_1_8 - 1.0) + 
	5.0d-1*(z4_1_8 - 1.0)*(z4_1_8 - 1.0) + 5.0d-1*(x5_1_8 - 1.0)*(x5_1_8 - 1.0) + 
	5.0d-1*(y5_1_8 - 1.0)*(y5_1_8 - 1.0) + 5.0d-1*(z5_1_8 - 1.0)*(z5_1_8 - 1.0) + 
	5.0d-1*(x6_1_8 - 1.0)*(x6_1_8 - 1.0) + 5.0d-1*(y6_1_8 - 1.0)*(y6_1_8 - 1.0) + 
	5.0d-1*(z6_1_8 - 1.0)*(z6_1_8 - 1.0) + 5.0d-1*(x7_1_8 - 1.0)*(x7_1_8 - 1.0) + 
	5.0d-1*(y7_1_8 - 1.0)*(y7_1_8 - 1.0) + 5.0d-1*(z7_1_8 - 1.0)*(z7_1_8 - 1.0) + 
	5.0d-1*(x8_1_8 - 1.0)*(x8_1_8 - 1.0) + 5.0d-1*(y8_1_8 - 1.0)*(y8_1_8 - 1.0) + 
	5.0d-1*(z8_1_8 - 1.0)*(z8_1_8 - 1.0) + 5.0d-1*(x9_1_8 - 1.0)*(x9_1_8 - 1.0) + 
	5.0d-1*(y9_1_8 - 1.0)*(y9_1_8 - 1.0) + 5.0d-1*(z9_1_8 - 1.0)*(z9_1_8 - 1.0) + 
	5.0d-1*(x1_2_8 - 1.0)*(x1_2_8 - 1.0) + 5.0d-1*(y1_2_8 - 1.0)*(y1_2_8 - 1.0) + 
	5.0d-1*(z1_2_8 - 1.0)*(z1_2_8 - 1.0) + 5.0d-1*(x2_2_8 - 1.0)*(x2_2_8 - 1.0) + 
	5.0d-1*(y2_2_8 - 1.0)*(y2_2_8 - 1.0) + 5.0d-1*(z2_2_8 - 1.0)*(z2_2_8 - 1.0) + 
	5.0d-1*(x3_2_8 - 1.0)*(x3_2_8 - 1.0) + 5.0d-1*(y3_2_8 - 1.0)*(y3_2_8 - 1.0) + 
	5.0d-1*(z3_2_8 - 1.0)*(z3_2_8 - 1.0) + 5.0d-1*(x4_2_8 - 1.0)*(x4_2_8 - 1.0) + 
	5.0d-1*(y4_2_8 - 1.0)*(y4_2_8 - 1.0) + 5.0d-1*(z4_2_8 - 1.0)*(z4_2_8 - 1.0) + 
	5.0d-1*(x5_2_8 - 1.0)*(x5_2_8 - 1.0) + 5.0d-1*(y5_2_8 - 1.0)*(y5_2_8 - 1.0) + 
	5.0d-1*(z5_2_8 - 1.0)*(z5_2_8 - 1.0) + 5.0d-1*(x6_2_8 - 1.0)*(x6_2_8 - 1.0) + 
	5.0d-1*(y6_2_8 - 1.0)*(y6_2_8 - 1.0) + 5.0d-1*(z6_2_8 - 1.0)*(z6_2_8 - 1.0) + 
	5.0d-1*(x7_2_8 - 1.0)*(x7_2_8 - 1.0) + 5.0d-1*(y7_2_8 - 1.0)*(y7_2_8 - 1.0) + 
	5.0d-1*(z7_2_8 - 1.0)*(z7_2_8 - 1.0) + 5.0d-1*(x8_2_8 - 1.0)*(x8_2_8 - 1.0) + 
	5.0d-1*(y8_2_8 - 1.0)*(y8_2_8 - 1.0) + 5.0d-1*(z8_2_8 - 1.0)*(z8_2_8 - 1.0) + 
	5.0d-1*(x9_2_8 - 1.0)*(x9_2_8 - 1.0) + 5.0d-1*(y9_2_8 - 1.0)*(y9_2_8 - 1.0) + 
	5.0d-1*(z9_2_8 - 1.0)*(z9_2_8 - 1.0) + 5.0d-1*(x1_3_8 - 1.0)*(x1_3_8 - 1.0) + 
	5.0d-1*(y1_3_8 - 1.0)*(y1_3_8 - 1.0) + 5.0d-1*(z1_3_8 - 1.0)*(z1_3_8 - 1.0) + 
	5.0d-1*(x2_3_8 - 1.0)*(x2_3_8 - 1.0) + 5.0d-1*(y2_3_8 - 1.0)*(y2_3_8 - 1.0) + 
	5.0d-1*(z2_3_8 - 1.0)*(z2_3_8 - 1.0) + 5.0d-1*(x3_3_8 - 1.0)*(x3_3_8 - 1.0) + 
	5.0d-1*(y3_3_8 - 1.0)*(y3_3_8 - 1.0) + 5.0d-1*(z3_3_8 - 1.0)*(z3_3_8 - 1.0) + 
	5.0d-1*(x4_3_8 - 1.0)*(x4_3_8 - 1.0) + 5.0d-1*(y4_3_8 - 1.0)*(y4_3_8 - 1.0) + 
	5.0d-1*(z4_3_8 - 1.0)*(z4_3_8 - 1.0) + 5.0d-1*(x5_3_8 - 1.0)*(x5_3_8 - 1.0) + 
	5.0d-1*(y5_3_8 - 1.0)*(y5_3_8 - 1.0) + 5.0d-1*(z5_3_8 - 1.0)*(z5_3_8 - 1.0) + 
	5.0d-1*(x6_3_8 - 1.0)*(x6_3_8 - 1.0) + 5.0d-1*(y6_3_8 - 1.0)*(y6_3_8 - 1.0) + 
	5.0d-1*(z6_3_8 - 1.0)*(z6_3_8 - 1.0) + 5.0d-1*(x7_3_8 - 1.0)*(x7_3_8 - 1.0) + 
	5.0d-1*(y7_3_8 - 1.0)*(y7_3_8 - 1.0) + 5.0d-1*(z7_3_8 - 1.0)*(z7_3_8 - 1.0) + 
	5.0d-1*(x8_3_8 - 1.0)*(x8_3_8 - 1.0) + 5.0d-1*(y8_3_8 - 1.0)*(y8_3_8 - 1.0) + 
	5.0d-1*(z8_3_8 - 1.0)*(z8_3_8 - 1.0) + 5.0d-1*(x9_3_8 - 1.0)*(x9_3_8 - 1.0) + 
	5.0d-1*(y9_3_8 - 1.0)*(y9_3_8 - 1.0) + 5.0d-1*(z9_3_8 - 1.0)*(z9_3_8 - 1.0) + 
	5.0d-1*(x1_4_8 - 1.0)*(x1_4_8 - 1.0) + 5.0d-1*(y1_4_8 - 1.0)*(y1_4_8 - 1.0) + 
	5.0d-1*(z1_4_8 - 1.0)*(z1_4_8 - 1.0) + 5.0d-1*(x2_4_8 - 1.0)*(x2_4_8 - 1.0) + 
	5.0d-1*(y2_4_8 - 1.0)*(y2_4_8 - 1.0) + 5.0d-1*(z2_4_8 - 1.0)*(z2_4_8 - 1.0) + 
	5.0d-1*(x3_4_8 - 1.0)*(x3_4_8 - 1.0) + 5.0d-1*(y3_4_8 - 1.0)*(y3_4_8 - 1.0) + 
	5.0d-1*(z3_4_8 - 1.0)*(z3_4_8 - 1.0) + 5.0d-1*(x4_4_8 - 1.0)*(x4_4_8 - 1.0) + 
	5.0d-1*(y4_4_8 - 1.0)*(y4_4_8 - 1.0) + 5.0d-1*(z4_4_8 - 1.0)*(z4_4_8 - 1.0) + 
	5.0d-1*(x5_4_8 - 1.0)*(x5_4_8 - 1.0) + 5.0d-1*(y5_4_8 - 1.0)*(y5_4_8 - 1.0) + 
	5.0d-1*(z5_4_8 - 1.0)*(z5_4_8 - 1.0) + 5.0d-1*(x6_4_8 - 1.0)*(x6_4_8 - 1.0) + 
	5.0d-1*(y6_4_8 - 1.0)*(y6_4_8 - 1.0) + 5.0d-1*(z6_4_8 - 1.0)*(z6_4_8 - 1.0) + 
	5.0d-1*(x7_4_8 - 1.0)*(x7_4_8 - 1.0) + 5.0d-1*(y7_4_8 - 1.0)*(y7_4_8 - 1.0) + 
	5.0d-1*(z7_4_8 - 1.0)*(z7_4_8 - 1.0) + 5.0d-1*(x8_4_8 - 1.0)*(x8_4_8 - 1.0) + 
	5.0d-1*(y8_4_8 - 1.0)*(y8_4_8 - 1.0) + 5.0d-1*(z8_4_8 - 1.0)*(z8_4_8 - 1.0) + 
	5.0d-1*(x9_4_8 - 1.0)*(x9_4_8 - 1.0) + 5.0d-1*(y9_4_8 - 1.0)*(y9_4_8 - 1.0) + 
	5.0d-1*(z9_4_8 - 1.0)*(z9_4_8 - 1.0) + 5.0d-1*(x1_5_8 - 1.0)*(x1_5_8 - 1.0) + 
	5.0d-1*(y1_5_8 - 1.0)*(y1_5_8 - 1.0) + 5.0d-1*(z1_5_8 - 1.0)*(z1_5_8 - 1.0) + 
	5.0d-1*(x2_5_8 - 1.0)*(x2_5_8 - 1.0) + 5.0d-1*(y2_5_8 - 1.0)*(y2_5_8 - 1.0) + 
	5.0d-1*(z2_5_8 - 1.0)*(z2_5_8 - 1.0) + 5.0d-1*(x3_5_8 - 1.0)*(x3_5_8 - 1.0) + 
	5.0d-1*(y3_5_8 - 1.0)*(y3_5_8 - 1.0) + 5.0d-1*(z3_5_8 - 1.0)*(z3_5_8 - 1.0) + 
	5.0d-1*(x4_5_8 - 1.0)*(x4_5_8 - 1.0) + 5.0d-1*(y4_5_8 - 1.0)*(y4_5_8 - 1.0) + 
	5.0d-1*(z4_5_8 - 1.0)*(z4_5_8 - 1.0) + 5.0d-1*(x5_5_8 - 1.0)*(x5_5_8 - 1.0) + 
	5.0d-1*(y5_5_8 - 1.0)*(y5_5_8 - 1.0) + 5.0d-1*(z5_5_8 - 1.0)*(z5_5_8 - 1.0) + 
	5.0d-1*(x6_5_8 - 1.0)*(x6_5_8 - 1.0) + 5.0d-1*(y6_5_8 - 1.0)*(y6_5_8 - 1.0) + 
	5.0d-1*(z6_5_8 - 1.0)*(z6_5_8 - 1.0) + 5.0d-1*(x7_5_8 - 1.0)*(x7_5_8 - 1.0) + 
	5.0d-1*(y7_5_8 - 1.0)*(y7_5_8 - 1.0) + 5.0d-1*(z7_5_8 - 1.0)*(z7_5_8 - 1.0) + 
	5.0d-1*(x8_5_8 - 1.0)*(x8_5_8 - 1.0) + 5.0d-1*(y8_5_8 - 1.0)*(y8_5_8 - 1.0) + 
	5.0d-1*(z8_5_8 - 1.0)*(z8_5_8 - 1.0) + 5.0d-1*(x9_5_8 - 1.0)*(x9_5_8 - 1.0) + 
	5.0d-1*(y9_5_8 - 1.0)*(y9_5_8 - 1.0) + 5.0d-1*(z9_5_8 - 1.0)*(z9_5_8 - 1.0) + 
	5.0d-1*(x1_6_8 - 1.0)*(x1_6_8 - 1.0) + 5.0d-1*(y1_6_8 - 1.0)*(y1_6_8 - 1.0) + 
	5.0d-1*(z1_6_8 - 1.0)*(z1_6_8 - 1.0) + 5.0d-1*(x2_6_8 - 1.0)*(x2_6_8 - 1.0) + 
	5.0d-1*(y2_6_8 - 1.0)*(y2_6_8 - 1.0) + 5.0d-1*(z2_6_8 - 1.0)*(z2_6_8 - 1.0) + 
	5.0d-1*(x3_6_8 - 1.0)*(x3_6_8 - 1.0) + 5.0d-1*(y3_6_8 - 1.0)*(y3_6_8 - 1.0) + 
	5.0d-1*(z3_6_8 - 1.0)*(z3_6_8 - 1.0) + 5.0d-1*(x4_6_8 - 1.0)*(x4_6_8 - 1.0) + 
	5.0d-1*(y4_6_8 - 1.0)*(y4_6_8 - 1.0) + 5.0d-1*(z4_6_8 - 1.0)*(z4_6_8 - 1.0) + 
	5.0d-1*(x5_6_8 - 1.0)*(x5_6_8 - 1.0) + 5.0d-1*(y5_6_8 - 1.0)*(y5_6_8 - 1.0) + 
	5.0d-1*(z5_6_8 - 1.0)*(z5_6_8 - 1.0) + 5.0d-1*(x6_6_8 - 1.0)*(x6_6_8 - 1.0) + 
	5.0d-1*(y6_6_8 - 1.0)*(y6_6_8 - 1.0) + 5.0d-1*(z6_6_8 - 1.0)*(z6_6_8 - 1.0) + 
	5.0d-1*(x7_6_8 - 1.0)*(x7_6_8 - 1.0) + 5.0d-1*(y7_6_8 - 1.0)*(y7_6_8 - 1.0) + 
	5.0d-1*(z7_6_8 - 1.0)*(z7_6_8 - 1.0) + 5.0d-1*(x8_6_8 - 1.0)*(x8_6_8 - 1.0) + 
	5.0d-1*(y8_6_8 - 1.0)*(y8_6_8 - 1.0) + 5.0d-1*(z8_6_8 - 1.0)*(z8_6_8 - 1.0) + 
	5.0d-1*(x9_6_8 - 1.0)*(x9_6_8 - 1.0) + 5.0d-1*(y9_6_8 - 1.0)*(y9_6_8 - 1.0) + 
	5.0d-1*(z9_6_8 - 1.0)*(z9_6_8 - 1.0) + 5.0d-1*(x1_7_8 - 1.0)*(x1_7_8 - 1.0) + 
	5.0d-1*(y1_7_8 - 1.0)*(y1_7_8 - 1.0) + 5.0d-1*(z1_7_8 - 1.0)*(z1_7_8 - 1.0) + 
	5.0d-1*(x2_7_8 - 1.0)*(x2_7_8 - 1.0) + 5.0d-1*(y2_7_8 - 1.0)*(y2_7_8 - 1.0) + 
	5.0d-1*(z2_7_8 - 1.0)*(z2_7_8 - 1.0) + 5.0d-1*(x3_7_8 - 1.0)*(x3_7_8 - 1.0) + 
	5.0d-1*(y3_7_8 - 1.0)*(y3_7_8 - 1.0) + 5.0d-1*(z3_7_8 - 1.0)*(z3_7_8 - 1.0) + 
	5.0d-1*(x4_7_8 - 1.0)*(x4_7_8 - 1.0) + 5.0d-1*(y4_7_8 - 1.0)*(y4_7_8 - 1.0) + 
	5.0d-1*(z4_7_8 - 1.0)*(z4_7_8 - 1.0) + 5.0d-1*(x5_7_8 - 1.0)*(x5_7_8 - 1.0) + 
	5.0d-1*(y5_7_8 - 1.0)*(y5_7_8 - 1.0) + 5.0d-1*(z5_7_8 - 1.0)*(z5_7_8 - 1.0) + 
	5.0d-1*(x6_7_8 - 1.0)*(x6_7_8 - 1.0) + 5.0d-1*(y6_7_8 - 1.0)*(y6_7_8 - 1.0) + 
	5.0d-1*(z6_7_8 - 1.0)*(z6_7_8 - 1.0) + 5.0d-1*(x7_7_8 - 1.0)*(x7_7_8 - 1.0) + 
	5.0d-1*(y7_7_8 - 1.0)*(y7_7_8 - 1.0) + 5.0d-1*(z7_7_8 - 1.0)*(z7_7_8 - 1.0) + 
	5.0d-1*(x8_7_8 - 1.0)*(x8_7_8 - 1.0) + 5.0d-1*(y8_7_8 - 1.0)*(y8_7_8 - 1.0) + 
	5.0d-1*(z8_7_8 - 1.0)*(z8_7_8 - 1.0) + 5.0d-1*(x9_7_8 - 1.0)*(x9_7_8 - 1.0) + 
	5.0d-1*(y9_7_8 - 1.0)*(y9_7_8 - 1.0) + 5.0d-1*(z9_7_8 - 1.0)*(z9_7_8 - 1.0) + 
	5.0d-1*(x1_8_8 - 1.0)*(x1_8_8 - 1.0) + 5.0d-1*(y1_8_8 - 1.0)*(y1_8_8 - 1.0) + 
	5.0d-1*(z1_8_8 - 1.0)*(z1_8_8 - 1.0) + 5.0d-1*(x2_8_8 - 1.0)*(x2_8_8 - 1.0) + 
	5.0d-1*(y2_8_8 - 1.0)*(y2_8_8 - 1.0) + 5.0d-1*(z2_8_8 - 1.0)*(z2_8_8 - 1.0) + 
	5.0d-1*(x3_8_8 - 1.0)*(x3_8_8 - 1.0) + 5.0d-1*(y3_8_8 - 1.0)*(y3_8_8 - 1.0) + 
	5.0d-1*(z3_8_8 - 1.0)*(z3_8_8 - 1.0) + 5.0d-1*(x4_8_8 - 1.0)*(x4_8_8 - 1.0) + 
	5.0d-1*(y4_8_8 - 1.0)*(y4_8_8 - 1.0) + 5.0d-1*(z4_8_8 - 1.0)*(z4_8_8 - 1.0) + 
	5.0d-1*(x5_8_8 - 1.0)*(x5_8_8 - 1.0) + 5.0d-1*(y5_8_8 - 1.0)*(y5_8_8 - 1.0) + 
	5.0d-1*(z5_8_8 - 1.0)*(z5_8_8 - 1.0) + 5.0d-1*(x6_8_8 - 1.0)*(x6_8_8 - 1.0) + 
	5.0d-1*(y6_8_8 - 1.0)*(y6_8_8 - 1.0) + 5.0d-1*(z6_8_8 - 1.0)*(z6_8_8 - 1.0) + 
	5.0d-1*(x7_8_8 - 1.0)*(x7_8_8 - 1.0) + 5.0d-1*(y7_8_8 - 1.0)*(y7_8_8 - 1.0) + 
	5.0d-1*(z7_8_8 - 1.0)*(z7_8_8 - 1.0) + 5.0d-1*(x8_8_8 - 1.0)*(x8_8_8 - 1.0) + 
	5.0d-1*(y8_8_8 - 1.0)*(y8_8_8 - 1.0) + 5.0d-1*(z8_8_8 - 1.0)*(z8_8_8 - 1.0) + 
	5.0d-1*(x9_8_8 - 1.0)*(x9_8_8 - 1.0) + 5.0d-1*(y9_8_8 - 1.0)*(y9_8_8 - 1.0) + 
	5.0d-1*(z9_8_8 - 1.0)*(z9_8_8 - 1.0) + 5.0d-1*(x1_9_8 - 1.0)*(x1_9_8 - 1.0) + 
	5.0d-1*(y1_9_8 - 1.0)*(y1_9_8 - 1.0) + 5.0d-1*(z1_9_8 - 1.0)*(z1_9_8 - 1.0) + 
	5.0d-1*(x2_9_8 - 1.0)*(x2_9_8 - 1.0) + 5.0d-1*(y2_9_8 - 1.0)*(y2_9_8 - 1.0) + 
	5.0d-1*(z2_9_8 - 1.0)*(z2_9_8 - 1.0) + 5.0d-1*(x3_9_8 - 1.0)*(x3_9_8 - 1.0) + 
	5.0d-1*(y3_9_8 - 1.0)*(y3_9_8 - 1.0) + 5.0d-1*(z3_9_8 - 1.0)*(z3_9_8 - 1.0) + 
	5.0d-1*(x4_9_8 - 1.0)*(x4_9_8 - 1.0) + 5.0d-1*(y4_9_8 - 1.0)*(y4_9_8 - 1.0) + 
	5.0d-1*(z4_9_8 - 1.0)*(z4_9_8 - 1.0) + 5.0d-1*(x5_9_8 - 1.0)*(x5_9_8 - 1.0) + 
	5.0d-1*(y5_9_8 - 1.0)*(y5_9_8 - 1.0) + 5.0d-1*(z5_9_8 - 1.0)*(z5_9_8 - 1.0) + 
	5.0d-1*(x6_9_8 - 1.0)*(x6_9_8 - 1.0) + 5.0d-1*(y6_9_8 - 1.0)*(y6_9_8 - 1.0) + 
	5.0d-1*(z6_9_8 - 1.0)*(z6_9_8 - 1.0) + 5.0d-1*(x7_9_8 - 1.0)*(x7_9_8 - 1.0) + 
	5.0d-1*(y7_9_8 - 1.0)*(y7_9_8 - 1.0) + 5.0d-1*(z7_9_8 - 1.0)*(z7_9_8 - 1.0) + 
	5.0d-1*(x8_9_8 - 1.0)*(x8_9_8 - 1.0) + 5.0d-1*(y8_9_8 - 1.0)*(y8_9_8 - 1.0) + 
	5.0d-1*(z8_9_8 - 1.0)*(z8_9_8 - 1.0) + 5.0d-1*(x9_9_8 - 1.0)*(x9_9_8 - 1.0) + 
	5.0d-1*(y9_9_8 - 1.0)*(y9_9_8 - 1.0) + 5.0d-1*(z9_9_8 - 1.0)*(z9_9_8 - 1.0) + 
	5.0d-1*(x1_1_9 - 1.0)*(x1_1_9 - 1.0) + 5.0d-1*(y1_1_9 - 1.0)*(y1_1_9 - 1.0) + 
	5.0d-1*(z1_1_9 - 1.0)*(z1_1_9 - 1.0) + 5.0d-1*(x2_1_9 - 1.0)*(x2_1_9 - 1.0) + 
	5.0d-1*(y2_1_9 - 1.0)*(y2_1_9 - 1.0) + 5.0d-1*(z2_1_9 - 1.0)*(z2_1_9 - 1.0) + 
	5.0d-1*(x3_1_9 - 1.0)*(x3_1_9 - 1.0) + 5.0d-1*(y3_1_9 - 1.0)*(y3_1_9 - 1.0) + 
	5.0d-1*(z3_1_9 - 1.0)*(z3_1_9 - 1.0) + 5.0d-1*(x4_1_9 - 1.0)*(x4_1_9 - 1.0) + 
	5.0d-1*(y4_1_9 - 1.0)*(y4_1_9 - 1.0) + 5.0d-1*(z4_1_9 - 1.0)*(z4_1_9 - 1.0) + 
	5.0d-1*(x5_1_9 - 1.0)*(x5_1_9 - 1.0) + 5.0d-1*(y5_1_9 - 1.0)*(y5_1_9 - 1.0) + 
	5.0d-1*(z5_1_9 - 1.0)*(z5_1_9 - 1.0) + 5.0d-1*(x6_1_9 - 1.0)*(x6_1_9 - 1.0) + 
	5.0d-1*(y6_1_9 - 1.0)*(y6_1_9 - 1.0) + 5.0d-1*(z6_1_9 - 1.0)*(z6_1_9 - 1.0) + 
	5.0d-1*(x7_1_9 - 1.0)*(x7_1_9 - 1.0) + 5.0d-1*(y7_1_9 - 1.0)*(y7_1_9 - 1.0) + 
	5.0d-1*(z7_1_9 - 1.0)*(z7_1_9 - 1.0) + 5.0d-1*(x8_1_9 - 1.0)*(x8_1_9 - 1.0) + 
	5.0d-1*(y8_1_9 - 1.0)*(y8_1_9 - 1.0) + 5.0d-1*(z8_1_9 - 1.0)*(z8_1_9 - 1.0) + 
	5.0d-1*(x9_1_9 - 1.0)*(x9_1_9 - 1.0) + 5.0d-1*(y9_1_9 - 1.0)*(y9_1_9 - 1.0) + 
	5.0d-1*(z9_1_9 - 1.0)*(z9_1_9 - 1.0) + 5.0d-1*(x1_2_9 - 1.0)*(x1_2_9 - 1.0) + 
	5.0d-1*(y1_2_9 - 1.0)*(y1_2_9 - 1.0) + 5.0d-1*(z1_2_9 - 1.0)*(z1_2_9 - 1.0) + 
	5.0d-1*(x2_2_9 - 1.0)*(x2_2_9 - 1.0) + 5.0d-1*(y2_2_9 - 1.0)*(y2_2_9 - 1.0) + 
	5.0d-1*(z2_2_9 - 1.0)*(z2_2_9 - 1.0) + 5.0d-1*(x3_2_9 - 1.0)*(x3_2_9 - 1.0) + 
	5.0d-1*(y3_2_9 - 1.0)*(y3_2_9 - 1.0) + 5.0d-1*(z3_2_9 - 1.0)*(z3_2_9 - 1.0) + 
	5.0d-1*(x4_2_9 - 1.0)*(x4_2_9 - 1.0) + 5.0d-1*(y4_2_9 - 1.0)*(y4_2_9 - 1.0) + 
	5.0d-1*(z4_2_9 - 1.0)*(z4_2_9 - 1.0) + 5.0d-1*(x5_2_9 - 1.0)*(x5_2_9 - 1.0) + 
	5.0d-1*(y5_2_9 - 1.0)*(y5_2_9 - 1.0) + 5.0d-1*(z5_2_9 - 1.0)*(z5_2_9 - 1.0) + 
	5.0d-1*(x6_2_9 - 1.0)*(x6_2_9 - 1.0) + 5.0d-1*(y6_2_9 - 1.0)*(y6_2_9 - 1.0) + 
	5.0d-1*(z6_2_9 - 1.0)*(z6_2_9 - 1.0) + 5.0d-1*(x7_2_9 - 1.0)*(x7_2_9 - 1.0) + 
	5.0d-1*(y7_2_9 - 1.0)*(y7_2_9 - 1.0) + 5.0d-1*(z7_2_9 - 1.0)*(z7_2_9 - 1.0) + 
	5.0d-1*(x8_2_9 - 1.0)*(x8_2_9 - 1.0) + 5.0d-1*(y8_2_9 - 1.0)*(y8_2_9 - 1.0) + 
	5.0d-1*(z8_2_9 - 1.0)*(z8_2_9 - 1.0) + 5.0d-1*(x9_2_9 - 1.0)*(x9_2_9 - 1.0) + 
	5.0d-1*(y9_2_9 - 1.0)*(y9_2_9 - 1.0) + 5.0d-1*(z9_2_9 - 1.0)*(z9_2_9 - 1.0) + 
	5.0d-1*(x1_3_9 - 1.0)*(x1_3_9 - 1.0) + 5.0d-1*(y1_3_9 - 1.0)*(y1_3_9 - 1.0) + 
	5.0d-1*(z1_3_9 - 1.0)*(z1_3_9 - 1.0) + 5.0d-1*(x2_3_9 - 1.0)*(x2_3_9 - 1.0) + 
	5.0d-1*(y2_3_9 - 1.0)*(y2_3_9 - 1.0) + 5.0d-1*(z2_3_9 - 1.0)*(z2_3_9 - 1.0) + 
	5.0d-1*(x3_3_9 - 1.0)*(x3_3_9 - 1.0) + 5.0d-1*(y3_3_9 - 1.0)*(y3_3_9 - 1.0) + 
	5.0d-1*(z3_3_9 - 1.0)*(z3_3_9 - 1.0) + 5.0d-1*(x4_3_9 - 1.0)*(x4_3_9 - 1.0) + 
	5.0d-1*(y4_3_9 - 1.0)*(y4_3_9 - 1.0) + 5.0d-1*(z4_3_9 - 1.0)*(z4_3_9 - 1.0) + 
	5.0d-1*(x5_3_9 - 1.0)*(x5_3_9 - 1.0) + 5.0d-1*(y5_3_9 - 1.0)*(y5_3_9 - 1.0) + 
	5.0d-1*(z5_3_9 - 1.0)*(z5_3_9 - 1.0) + 5.0d-1*(x6_3_9 - 1.0)*(x6_3_9 - 1.0) + 
	5.0d-1*(y6_3_9 - 1.0)*(y6_3_9 - 1.0) + 5.0d-1*(z6_3_9 - 1.0)*(z6_3_9 - 1.0) + 
	5.0d-1*(x7_3_9 - 1.0)*(x7_3_9 - 1.0) + 5.0d-1*(y7_3_9 - 1.0)*(y7_3_9 - 1.0) + 
	5.0d-1*(z7_3_9 - 1.0)*(z7_3_9 - 1.0) + 5.0d-1*(x8_3_9 - 1.0)*(x8_3_9 - 1.0) + 
	5.0d-1*(y8_3_9 - 1.0)*(y8_3_9 - 1.0) + 5.0d-1*(z8_3_9 - 1.0)*(z8_3_9 - 1.0) + 
	5.0d-1*(x9_3_9 - 1.0)*(x9_3_9 - 1.0) + 5.0d-1*(y9_3_9 - 1.0)*(y9_3_9 - 1.0) + 
	5.0d-1*(z9_3_9 - 1.0)*(z9_3_9 - 1.0) + 5.0d-1*(x1_4_9 - 1.0)*(x1_4_9 - 1.0) + 
	5.0d-1*(y1_4_9 - 1.0)*(y1_4_9 - 1.0) + 5.0d-1*(z1_4_9 - 1.0)*(z1_4_9 - 1.0) + 
	5.0d-1*(x2_4_9 - 1.0)*(x2_4_9 - 1.0) + 5.0d-1*(y2_4_9 - 1.0)*(y2_4_9 - 1.0) + 
	5.0d-1*(z2_4_9 - 1.0)*(z2_4_9 - 1.0) + 5.0d-1*(x3_4_9 - 1.0)*(x3_4_9 - 1.0) + 
	5.0d-1*(y3_4_9 - 1.0)*(y3_4_9 - 1.0) + 5.0d-1*(z3_4_9 - 1.0)*(z3_4_9 - 1.0) + 
	5.0d-1*(x4_4_9 - 1.0)*(x4_4_9 - 1.0) + 5.0d-1*(y4_4_9 - 1.0)*(y4_4_9 - 1.0) + 
	5.0d-1*(z4_4_9 - 1.0)*(z4_4_9 - 1.0) + 5.0d-1*(x5_4_9 - 1.0)*(x5_4_9 - 1.0) + 
	5.0d-1*(y5_4_9 - 1.0)*(y5_4_9 - 1.0) + 5.0d-1*(z5_4_9 - 1.0)*(z5_4_9 - 1.0) + 
	5.0d-1*(x6_4_9 - 1.0)*(x6_4_9 - 1.0) + 5.0d-1*(y6_4_9 - 1.0)*(y6_4_9 - 1.0) + 
	5.0d-1*(z6_4_9 - 1.0)*(z6_4_9 - 1.0) + 5.0d-1*(x7_4_9 - 1.0)*(x7_4_9 - 1.0) + 
	5.0d-1*(y7_4_9 - 1.0)*(y7_4_9 - 1.0) + 5.0d-1*(z7_4_9 - 1.0)*(z7_4_9 - 1.0) + 
	5.0d-1*(x8_4_9 - 1.0)*(x8_4_9 - 1.0) + 5.0d-1*(y8_4_9 - 1.0)*(y8_4_9 - 1.0) + 
	5.0d-1*(z8_4_9 - 1.0)*(z8_4_9 - 1.0) + 5.0d-1*(x9_4_9 - 1.0)*(x9_4_9 - 1.0) + 
	5.0d-1*(y9_4_9 - 1.0)*(y9_4_9 - 1.0) + 5.0d-1*(z9_4_9 - 1.0)*(z9_4_9 - 1.0) + 
	5.0d-1*(x1_5_9 - 1.0)*(x1_5_9 - 1.0) + 5.0d-1*(y1_5_9 - 1.0)*(y1_5_9 - 1.0) + 
	5.0d-1*(z1_5_9 - 1.0)*(z1_5_9 - 1.0) + 5.0d-1*(x2_5_9 - 1.0)*(x2_5_9 - 1.0) + 
	5.0d-1*(y2_5_9 - 1.0)*(y2_5_9 - 1.0) + 5.0d-1*(z2_5_9 - 1.0)*(z2_5_9 - 1.0) + 
	5.0d-1*(x3_5_9 - 1.0)*(x3_5_9 - 1.0) + 5.0d-1*(y3_5_9 - 1.0)*(y3_5_9 - 1.0) + 
	5.0d-1*(z3_5_9 - 1.0)*(z3_5_9 - 1.0) + 5.0d-1*(x4_5_9 - 1.0)*(x4_5_9 - 1.0) + 
	5.0d-1*(y4_5_9 - 1.0)*(y4_5_9 - 1.0) + 5.0d-1*(z4_5_9 - 1.0)*(z4_5_9 - 1.0) + 
	5.0d-1*(x5_5_9 - 1.0)*(x5_5_9 - 1.0) + 5.0d-1*(y5_5_9 - 1.0)*(y5_5_9 - 1.0) + 
	5.0d-1*(z5_5_9 - 1.0)*(z5_5_9 - 1.0) + 5.0d-1*(x6_5_9 - 1.0)*(x6_5_9 - 1.0) + 
	5.0d-1*(y6_5_9 - 1.0)*(y6_5_9 - 1.0) + 5.0d-1*(z6_5_9 - 1.0)*(z6_5_9 - 1.0) + 
	5.0d-1*(x7_5_9 - 1.0)*(x7_5_9 - 1.0) + 5.0d-1*(y7_5_9 - 1.0)*(y7_5_9 - 1.0) + 
	5.0d-1*(z7_5_9 - 1.0)*(z7_5_9 - 1.0) + 5.0d-1*(x8_5_9 - 1.0)*(x8_5_9 - 1.0) + 
	5.0d-1*(y8_5_9 - 1.0)*(y8_5_9 - 1.0) + 5.0d-1*(z8_5_9 - 1.0)*(z8_5_9 - 1.0) + 
	5.0d-1*(x9_5_9 - 1.0)*(x9_5_9 - 1.0) + 5.0d-1*(y9_5_9 - 1.0)*(y9_5_9 - 1.0) + 
	5.0d-1*(z9_5_9 - 1.0)*(z9_5_9 - 1.0) + 5.0d-1*(x1_6_9 - 1.0)*(x1_6_9 - 1.0) + 
	5.0d-1*(y1_6_9 - 1.0)*(y1_6_9 - 1.0) + 5.0d-1*(z1_6_9 - 1.0)*(z1_6_9 - 1.0) + 
	5.0d-1*(x2_6_9 - 1.0)*(x2_6_9 - 1.0) + 5.0d-1*(y2_6_9 - 1.0)*(y2_6_9 - 1.0) + 
	5.0d-1*(z2_6_9 - 1.0)*(z2_6_9 - 1.0) + 5.0d-1*(x3_6_9 - 1.0)*(x3_6_9 - 1.0) + 
	5.0d-1*(y3_6_9 - 1.0)*(y3_6_9 - 1.0) + 5.0d-1*(z3_6_9 - 1.0)*(z3_6_9 - 1.0) + 
	5.0d-1*(x4_6_9 - 1.0)*(x4_6_9 - 1.0) + 5.0d-1*(y4_6_9 - 1.0)*(y4_6_9 - 1.0) + 
	5.0d-1*(z4_6_9 - 1.0)*(z4_6_9 - 1.0) + 5.0d-1*(x5_6_9 - 1.0)*(x5_6_9 - 1.0) + 
	5.0d-1*(y5_6_9 - 1.0)*(y5_6_9 - 1.0) + 5.0d-1*(z5_6_9 - 1.0)*(z5_6_9 - 1.0) + 
	5.0d-1*(x6_6_9 - 1.0)*(x6_6_9 - 1.0) + 5.0d-1*(y6_6_9 - 1.0)*(y6_6_9 - 1.0) + 
	5.0d-1*(z6_6_9 - 1.0)*(z6_6_9 - 1.0) + 5.0d-1*(x7_6_9 - 1.0)*(x7_6_9 - 1.0) + 
	5.0d-1*(y7_6_9 - 1.0)*(y7_6_9 - 1.0) + 5.0d-1*(z7_6_9 - 1.0)*(z7_6_9 - 1.0) + 
	5.0d-1*(x8_6_9 - 1.0)*(x8_6_9 - 1.0) + 5.0d-1*(y8_6_9 - 1.0)*(y8_6_9 - 1.0) + 
	5.0d-1*(z8_6_9 - 1.0)*(z8_6_9 - 1.0) + 5.0d-1*(x9_6_9 - 1.0)*(x9_6_9 - 1.0) + 
	5.0d-1*(y9_6_9 - 1.0)*(y9_6_9 - 1.0) + 5.0d-1*(z9_6_9 - 1.0)*(z9_6_9 - 1.0) + 
	5.0d-1*(x1_7_9 - 1.0)*(x1_7_9 - 1.0) + 5.0d-1*(y1_7_9 - 1.0)*(y1_7_9 - 1.0) + 
	5.0d-1*(z1_7_9 - 1.0)*(z1_7_9 - 1.0) + 5.0d-1*(x2_7_9 - 1.0)*(x2_7_9 - 1.0) + 
	5.0d-1*(y2_7_9 - 1.0)*(y2_7_9 - 1.0) + 5.0d-1*(z2_7_9 - 1.0)*(z2_7_9 - 1.0) + 
	5.0d-1*(x3_7_9 - 1.0)*(x3_7_9 - 1.0) + 5.0d-1*(y3_7_9 - 1.0)*(y3_7_9 - 1.0) + 
	5.0d-1*(z3_7_9 - 1.0)*(z3_7_9 - 1.0) + 5.0d-1*(x4_7_9 - 1.0)*(x4_7_9 - 1.0) + 
	5.0d-1*(y4_7_9 - 1.0)*(y4_7_9 - 1.0) + 5.0d-1*(z4_7_9 - 1.0)*(z4_7_9 - 1.0) + 
	5.0d-1*(x5_7_9 - 1.0)*(x5_7_9 - 1.0) + 5.0d-1*(y5_7_9 - 1.0)*(y5_7_9 - 1.0) + 
	5.0d-1*(z5_7_9 - 1.0)*(z5_7_9 - 1.0) + 5.0d-1*(x6_7_9 - 1.0)*(x6_7_9 - 1.0) + 
	5.0d-1*(y6_7_9 - 1.0)*(y6_7_9 - 1.0) + 5.0d-1*(z6_7_9 - 1.0)*(z6_7_9 - 1.0) + 
	5.0d-1*(x7_7_9 - 1.0)*(x7_7_9 - 1.0) + 5.0d-1*(y7_7_9 - 1.0)*(y7_7_9 - 1.0) + 
	5.0d-1*(z7_7_9 - 1.0)*(z7_7_9 - 1.0) + 5.0d-1*(x8_7_9 - 1.0)*(x8_7_9 - 1.0) + 
	5.0d-1*(y8_7_9 - 1.0)*(y8_7_9 - 1.0) + 5.0d-1*(z8_7_9 - 1.0)*(z8_7_9 - 1.0) + 
	5.0d-1*(x9_7_9 - 1.0)*(x9_7_9 - 1.0) + 5.0d-1*(y9_7_9 - 1.0)*(y9_7_9 - 1.0) + 
	5.0d-1*(z9_7_9 - 1.0)*(z9_7_9 - 1.0) + 5.0d-1*(x1_8_9 - 1.0)*(x1_8_9 - 1.0) + 
	5.0d-1*(y1_8_9 - 1.0)*(y1_8_9 - 1.0) + 5.0d-1*(z1_8_9 - 1.0)*(z1_8_9 - 1.0) + 
	5.0d-1*(x2_8_9 - 1.0)*(x2_8_9 - 1.0) + 5.0d-1*(y2_8_9 - 1.0)*(y2_8_9 - 1.0) + 
	5.0d-1*(z2_8_9 - 1.0)*(z2_8_9 - 1.0) + 5.0d-1*(x3_8_9 - 1.0)*(x3_8_9 - 1.0) + 
	5.0d-1*(y3_8_9 - 1.0)*(y3_8_9 - 1.0) + 5.0d-1*(z3_8_9 - 1.0)*(z3_8_9 - 1.0) + 
	5.0d-1*(x4_8_9 - 1.0)*(x4_8_9 - 1.0) + 5.0d-1*(y4_8_9 - 1.0)*(y4_8_9 - 1.0) + 
	5.0d-1*(z4_8_9 - 1.0)*(z4_8_9 - 1.0) + 5.0d-1*(x5_8_9 - 1.0)*(x5_8_9 - 1.0) + 
	5.0d-1*(y5_8_9 - 1.0)*(y5_8_9 - 1.0) + 5.0d-1*(z5_8_9 - 1.0)*(z5_8_9 - 1.0) + 
	5.0d-1*(x6_8_9 - 1.0)*(x6_8_9 - 1.0) + 5.0d-1*(y6_8_9 - 1.0)*(y6_8_9 - 1.0) + 
	5.0d-1*(z6_8_9 - 1.0)*(z6_8_9 - 1.0) + 5.0d-1*(x7_8_9 - 1.0)*(x7_8_9 - 1.0) + 
	5.0d-1*(y7_8_9 - 1.0)*(y7_8_9 - 1.0) + 5.0d-1*(z7_8_9 - 1.0)*(z7_8_9 - 1.0) + 
	5.0d-1*(x8_8_9 - 1.0)*(x8_8_9 - 1.0) + 5.0d-1*(y8_8_9 - 1.0)*(y8_8_9 - 1.0) + 
	5.0d-1*(z8_8_9 - 1.0)*(z8_8_9 - 1.0) + 5.0d-1*(x9_8_9 - 1.0)*(x9_8_9 - 1.0) + 
	5.0d-1*(y9_8_9 - 1.0)*(y9_8_9 - 1.0) + 5.0d-1*(z9_8_9 - 1.0)*(z9_8_9 - 1.0) + 
	5.0d-1*(x1_9_9 - 1.0)*(x1_9_9 - 1.0) + 5.0d-1*(y1_9_9 - 1.0)*(y1_9_9 - 1.0) + 
	5.0d-1*(z1_9_9 - 1.0)*(z1_9_9 - 1.0) + 5.0d-1*(x2_9_9 - 1.0)*(x2_9_9 - 1.0) + 
	5.0d-1*(y2_9_9 - 1.0)*(y2_9_9 - 1.0) + 5.0d-1*(z2_9_9 - 1.0)*(z2_9_9 - 1.0) + 
	5.0d-1*(x3_9_9 - 1.0)*(x3_9_9 - 1.0) + 5.0d-1*(y3_9_9 - 1.0)*(y3_9_9 - 1.0) + 
	5.0d-1*(z3_9_9 - 1.0)*(z3_9_9 - 1.0) + 5.0d-1*(x4_9_9 - 1.0)*(x4_9_9 - 1.0) + 
	5.0d-1*(y4_9_9 - 1.0)*(y4_9_9 - 1.0) + 5.0d-1*(z4_9_9 - 1.0)*(z4_9_9 - 1.0) + 
	5.0d-1*(x5_9_9 - 1.0)*(x5_9_9 - 1.0) + 5.0d-1*(y5_9_9 - 1.0)*(y5_9_9 - 1.0) + 
	5.0d-1*(z5_9_9 - 1.0)*(z5_9_9 - 1.0) + 5.0d-1*(x6_9_9 - 1.0)*(x6_9_9 - 1.0) + 
	5.0d-1*(y6_9_9 - 1.0)*(y6_9_9 - 1.0) + 5.0d-1*(z6_9_9 - 1.0)*(z6_9_9 - 1.0) + 
	5.0d-1*(x7_9_9 - 1.0)*(x7_9_9 - 1.0) + 5.0d-1*(y7_9_9 - 1.0)*(y7_9_9 - 1.0) + 
	5.0d-1*(z7_9_9 - 1.0)*(z7_9_9 - 1.0) + 5.0d-1*(x8_9_9 - 1.0)*(x8_9_9 - 1.0) + 
	5.0d-1*(y8_9_9 - 1.0)*(y8_9_9 - 1.0) + 5.0d-1*(z8_9_9 - 1.0)*(z8_9_9 - 1.0) + 
	5.0d-1*(x9_9_9 - 1.0)*(x9_9_9 - 1.0) + 5.0d-1*(y9_9_9 - 1.0)*(y9_9_9 - 1.0) + 
	5.0d-1*(z9_9_9 - 1.0)*(z9_9_9 - 1.0) + 5.0d-1*(y10_1_1 - 1.0)*(y10_1_1 - 1.0) + 
	5.0d-1*(z10_1_1 - 1.0)*(z10_1_1 - 1.0) + 5.0d-1*(y10_2_1 - 1.0)*(y10_2_1 - 1.0) 
	+ 5.0d-1*(z10_2_1 - 1.0)*(z10_2_1 - 1.0) + 5.0d-1*(y10_3_1 - 1.0)*(y10_3_1 - 
	1.0) + 5.0d-1*(z10_3_1 - 1.0)*(z10_3_1 - 1.0) + 5.0d-1*(y10_4_1 - 1.0)*(y10_4_1 
	- 1.0) + 5.0d-1*(z10_4_1 - 1.0)*(z10_4_1 - 1.0) + 5.0d-1*(y10_5_1 - 
	1.0)*(y10_5_1 - 1.0) + 5.0d-1*(z10_5_1 - 1.0)*(z10_5_1 - 1.0) + 5.0d-1*(y10_6_1 
	- 1.0)*(y10_6_1 - 1.0) + 5.0d-1*(z10_6_1 - 1.0)*(z10_6_1 - 1.0) + 
	5.0d-1*(y10_7_1 - 1.0)*(y10_7_1 - 1.0) + 5.0d-1*(z10_7_1 - 1.0)*(z10_7_1 - 1.0) 
	+ 5.0d-1*(y10_8_1 - 1.0)*(y10_8_1 - 1.0) + 5.0d-1*(z10_8_1 - 1.0)*(z10_8_1 - 
	1.0) + 5.0d-1*(y10_9_1 - 1.0)*(y10_9_1 - 1.0) + 5.0d-1*(z10_9_1 - 1.0)*(z10_9_1 
	- 1.0) + 5.0d-1*(y10_1_2 - 1.0)*(y10_1_2 - 1.0) + 5.0d-1*(z10_1_2 - 
	1.0)*(z10_1_2 - 1.0) + 5.0d-1*(y10_2_2 - 1.0)*(y10_2_2 - 1.0) + 5.0d-1*(z10_2_2 
	- 1.0)*(z10_2_2 - 1.0) + 5.0d-1*(y10_3_2 - 1.0)*(y10_3_2 - 1.0) + 
	5.0d-1*(z10_3_2 - 1.0)*(z10_3_2 - 1.0) + 5.0d-1*(y10_4_2 - 1.0)*(y10_4_2 - 1.0) 
	+ 5.0d-1*(z10_4_2 - 1.0)*(z10_4_2 - 1.0) + 5.0d-1*(y10_5_2 - 1.0)*(y10_5_2 - 
	1.0) + 5.0d-1*(z10_5_2 - 1.0)*(z10_5_2 - 1.0) + 5.0d-1*(y10_6_2 - 1.0)*(y10_6_2 
	- 1.0) + 5.0d-1*(z10_6_2 - 1.0)*(z10_6_2 - 1.0) + 5.0d-1*(y10_7_2 - 
	1.0)*(y10_7_2 - 1.0) + 5.0d-1*(z10_7_2 - 1.0)*(z10_7_2 - 1.0) + 5.0d-1*(y10_8_2 
	- 1.0)*(y10_8_2 - 1.0) + 5.0d-1*(z10_8_2 - 1.0)*(z10_8_2 - 1.0) + 
	5.0d-1*(y10_9_2 - 1.0)*(y10_9_2 - 1.0) + 5.0d-1*(z10_9_2 - 1.0)*(z10_9_2 - 1.0) 
	+ 5.0d-1*(y10_1_3 - 1.0)*(y10_1_3 - 1.0) + 5.0d-1*(z10_1_3 - 1.0)*(z10_1_3 - 
	1.0) + 5.0d-1*(y10_2_3 - 1.0)*(y10_2_3 - 1.0) + 5.0d-1*(z10_2_3 - 1.0)*(z10_2_3 
	- 1.0) + 5.0d-1*(y10_3_3 - 1.0)*(y10_3_3 - 1.0) + 5.0d-1*(z10_3_3 - 
	1.0)*(z10_3_3 - 1.0) + 5.0d-1*(y10_4_3 - 1.0)*(y10_4_3 - 1.0) + 5.0d-1*(z10_4_3 
	- 1.0)*(z10_4_3 - 1.0) + 5.0d-1*(y10_5_3 - 1.0)*(y10_5_3 - 1.0) + 
	5.0d-1*(z10_5_3 - 1.0)*(z10_5_3 - 1.0) + 5.0d-1*(y10_6_3 - 1.0)*(y10_6_3 - 1.0) 
	+ 5.0d-1*(z10_6_3 - 1.0)*(z10_6_3 - 1.0) + 5.0d-1*(y10_7_3 - 1.0)*(y10_7_3 - 
	1.0) + 5.0d-1*(z10_7_3 - 1.0)*(z10_7_3 - 1.0) + 5.0d-1*(y10_8_3 - 1.0)*(y10_8_3 
	- 1.0) + 5.0d-1*(z10_8_3 - 1.0)*(z10_8_3 - 1.0) + 5.0d-1*(y10_9_3 - 
	1.0)*(y10_9_3 - 1.0) + 5.0d-1*(z10_9_3 - 1.0)*(z10_9_3 - 1.0) + 5.0d-1*(y10_1_4 
	- 1.0)*(y10_1_4 - 1.0) + 5.0d-1*(z10_1_4 - 1.0)*(z10_1_4 - 1.0) + 
	5.0d-1*(y10_2_4 - 1.0)*(y10_2_4 - 1.0) + 5.0d-1*(z10_2_4 - 1.0)*(z10_2_4 - 1.0) 
	+ 5.0d-1*(y10_3_4 - 1.0)*(y10_3_4 - 1.0) + 5.0d-1*(z10_3_4 - 1.0)*(z10_3_4 - 
	1.0) + 5.0d-1*(y10_4_4 - 1.0)*(y10_4_4 - 1.0) + 5.0d-1*(z10_4_4 - 1.0)*(z10_4_4 
	- 1.0) + 5.0d-1*(y10_5_4 - 1.0)*(y10_5_4 - 1.0) + 5.0d-1*(z10_5_4 - 
	1.0)*(z10_5_4 - 1.0) + 5.0d-1*(y10_6_4 - 1.0)*(y10_6_4 - 1.0) + 5.0d-1*(z10_6_4 
	- 1.0)*(z10_6_4 - 1.0) + 5.0d-1*(y10_7_4 - 1.0)*(y10_7_4 - 1.0) + 
	5.0d-1*(z10_7_4 - 1.0)*(z10_7_4 - 1.0) + 5.0d-1*(y10_8_4 - 1.0)*(y10_8_4 - 1.0) 
	+ 5.0d-1*(z10_8_4 - 1.0)*(z10_8_4 - 1.0) + 5.0d-1*(y10_9_4 - 1.0)*(y10_9_4 - 
	1.0) + 5.0d-1*(z10_9_4 - 1.0)*(z10_9_4 - 1.0) + 5.0d-1*(y10_1_5 - 1.0)*(y10_1_5 
	- 1.0) + 5.0d-1*(z10_1_5 - 1.0)*(z10_1_5 - 1.0) + 5.0d-1*(y10_2_5 - 
	1.0)*(y10_2_5 - 1.0) + 5.0d-1*(z10_2_5 - 1.0)*(z10_2_5 - 1.0) + 5.0d-1*(y10_3_5 
	- 1.0)*(y10_3_5 - 1.0) + 5.0d-1*(z10_3_5 - 1.0)*(z10_3_5 - 1.0) + 
	5.0d-1*(y10_4_5 - 1.0)*(y10_4_5 - 1.0) + 5.0d-1*(z10_4_5 - 1.0)*(z10_4_5 - 1.0) 
	+ 5.0d-1*(y10_5_5 - 1.0)*(y10_5_5 - 1.0) + 5.0d-1*(z10_5_5 - 1.0)*(z10_5_5 - 
	1.0) + 5.0d-1*(y10_6_5 - 1.0)*(y10_6_5 - 1.0) + 5.0d-1*(z10_6_5 - 1.0)*(z10_6_5 
	- 1.0) + 5.0d-1*(y10_7_5 - 1.0)*(y10_7_5 - 1.0) + 5.0d-1*(z10_7_5 - 
	1.0)*(z10_7_5 - 1.0) + 5.0d-1*(y10_8_5 - 1.0)*(y10_8_5 - 1.0) + 5.0d-1*(z10_8_5 
	- 1.0)*(z10_8_5 - 1.0) + 5.0d-1*(y10_9_5 - 1.0)*(y10_9_5 - 1.0) + 
	5.0d-1*(z10_9_5 - 1.0)*(z10_9_5 - 1.0) + 5.0d-1*(y10_1_6 - 1.0)*(y10_1_6 - 1.0) 
	+ 5.0d-1*(z10_1_6 - 1.0)*(z10_1_6 - 1.0) + 5.0d-1*(y10_2_6 - 1.0)*(y10_2_6 - 
	1.0) + 5.0d-1*(z10_2_6 - 1.0)*(z10_2_6 - 1.0) + 5.0d-1*(y10_3_6 - 1.0)*(y10_3_6 
	- 1.0) + 5.0d-1*(z10_3_6 - 1.0)*(z10_3_6 - 1.0) + 5.0d-1*(y10_4_6 - 
	1.0)*(y10_4_6 - 1.0) + 5.0d-1*(z10_4_6 - 1.0)*(z10_4_6 - 1.0) + 5.0d-1*(y10_5_6 
	- 1.0)*(y10_5_6 - 1.0) + 5.0d-1*(z10_5_6 - 1.0)*(z10_5_6 - 1.0) + 
	5.0d-1*(y10_6_6 - 1.0)*(y10_6_6 - 1.0) + 5.0d-1*(z10_6_6 - 1.0)*(z10_6_6 - 1.0) 
	+ 5.0d-1*(y10_7_6 - 1.0)*(y10_7_6 - 1.0) + 5.0d-1*(z10_7_6 - 1.0)*(z10_7_6 - 
	1.0) + 5.0d-1*(y10_8_6 - 1.0)*(y10_8_6 - 1.0) + 5.0d-1*(z10_8_6 - 1.0)*(z10_8_6 
	- 1.0) + 5.0d-1*(y10_9_6 - 1.0)*(y10_9_6 - 1.0) + 5.0d-1*(z10_9_6 - 
	1.0)*(z10_9_6 - 1.0) + 5.0d-1*(y10_1_7 - 1.0)*(y10_1_7 - 1.0) + 5.0d-1*(z10_1_7 
	- 1.0)*(z10_1_7 - 1.0) + 5.0d-1*(y10_2_7 - 1.0)*(y10_2_7 - 1.0) + 
	5.0d-1*(z10_2_7 - 1.0)*(z10_2_7 - 1.0) + 5.0d-1*(y10_3_7 - 1.0)*(y10_3_7 - 1.0) 
	+ 5.0d-1*(z10_3_7 - 1.0)*(z10_3_7 - 1.0) + 5.0d-1*(y10_4_7 - 1.0)*(y10_4_7 - 
	1.0) + 5.0d-1*(z10_4_7 - 1.0)*(z10_4_7 - 1.0) + 5.0d-1*(y10_5_7 - 1.0)*(y10_5_7 
	- 1.0) + 5.0d-1*(z10_5_7 - 1.0)*(z10_5_7 - 1.0) + 5.0d-1*(y10_6_7 - 
	1.0)*(y10_6_7 - 1.0) + 5.0d-1*(z10_6_7 - 1.0)*(z10_6_7 - 1.0) + 5.0d-1*(y10_7_7 
	- 1.0)*(y10_7_7 - 1.0) + 5.0d-1*(z10_7_7 - 1.0)*(z10_7_7 - 1.0) + 
	5.0d-1*(y10_8_7 - 1.0)*(y10_8_7 - 1.0) + 5.0d-1*(z10_8_7 - 1.0)*(z10_8_7 - 1.0) 
	+ 5.0d-1*(y10_9_7 - 1.0)*(y10_9_7 - 1.0) + 5.0d-1*(z10_9_7 - 1.0)*(z10_9_7 - 
	1.0) + 5.0d-1*(y10_1_8 - 1.0)*(y10_1_8 - 1.0) + 5.0d-1*(z10_1_8 - 1.0)*(z10_1_8 
	- 1.0) + 5.0d-1*(y10_2_8 - 1.0)*(y10_2_8 - 1.0) + 5.0d-1*(z10_2_8 - 
	1.0)*(z10_2_8 - 1.0) + 5.0d-1*(y10_3_8 - 1.0)*(y10_3_8 - 1.0) + 5.0d-1*(z10_3_8 
	- 1.0)*(z10_3_8 - 1.0) + 5.0d-1*(y10_4_8 - 1.0)*(y10_4_8 - 1.0) + 
	5.0d-1*(z10_4_8 - 1.0)*(z10_4_8 - 1.0) + 5.0d-1*(y10_5_8 - 1.0)*(y10_5_8 - 1.0) 
	+ 5.0d-1*(z10_5_8 - 1.0)*(z10_5_8 - 1.0) + 5.0d-1*(y10_6_8 - 1.0)*(y10_6_8 - 
	1.0) + 5.0d-1*(z10_6_8 - 1.0)*(z10_6_8 - 1.0) + 5.0d-1*(y10_7_8 - 1.0)*(y10_7_8 
	- 1.0) + 5.0d-1*(z10_7_8 - 1.0)*(z10_7_8 - 1.0) + 5.0d-1*(y10_8_8 - 
	1.0)*(y10_8_8 - 1.0) + 5.0d-1*(z10_8_8 - 1.0)*(z10_8_8 - 1.0) + 5.0d-1*(y10_9_8 
	- 1.0)*(y10_9_8 - 1.0) + 5.0d-1*(z10_9_8 - 1.0)*(z10_9_8 - 1.0) + 
	5.0d-1*(y10_1_9 - 1.0)*(y10_1_9 - 1.0) + 5.0d-1*(z10_1_9 - 1.0)*(z10_1_9 - 1.0) 
	+ 5.0d-1*(y10_2_9 - 1.0)*(y10_2_9 - 1.0) + 5.0d-1*(z10_2_9 - 1.0)*(z10_2_9 - 
	1.0) + 5.0d-1*(y10_3_9 - 1.0)*(y10_3_9 - 1.0) + 5.0d-1*(z10_3_9 - 1.0)*(z10_3_9 
	- 1.0) + 5.0d-1*(y10_4_9 - 1.0)*(y10_4_9 - 1.0) + 5.0d-1*(z10_4_9 - 
	1.0)*(z10_4_9 - 1.0) + 5.0d-1*(y10_5_9 - 1.0)*(y10_5_9 - 1.0) + 5.0d-1*(z10_5_9 
	- 1.0)*(z10_5_9 - 1.0) + 5.0d-1*(y10_6_9 - 1.0)*(y10_6_9 - 1.0) + 
	5.0d-1*(z10_6_9 - 1.0)*(z10_6_9 - 1.0) + 5.0d-1*(y10_7_9 - 1.0)*(y10_7_9 - 1.0) 
	+ 5.0d-1*(z10_7_9 - 1.0)*(z10_7_9 - 1.0) + 5.0d-1*(y10_8_9 - 1.0)*(y10_8_9 - 
	1.0) + 5.0d-1*(z10_8_9 - 1.0)*(z10_8_9 - 1.0) + 5.0d-1*(y10_9_9 - 1.0)*(y10_9_9 
	- 1.0) + 5.0d-1*(z10_9_9 - 1.0)*(z10_9_9 - 1.0) + 5.0d-1*(x1_10_1 - 
	1.0)*(x1_10_1 - 1.0) + 5.0d-1*(z1_10_1 - 1.0)*(z1_10_1 - 1.0) + 5.0d-1*(x2_10_1 
	- 1.0)*(x2_10_1 - 1.0) + 5.0d-1*(z2_10_1 - 1.0)*(z2_10_1 - 1.0) + 
	5.0d-1*(x3_10_1 - 1.0)*(x3_10_1 - 1.0) + 5.0d-1*(z3_10_1 - 1.0)*(z3_10_1 - 1.0) 
	+ 5.0d-1*(x4_10_1 - 1.0)*(x4_10_1 - 1.0) + 5.0d-1*(z4_10_1 - 1.0)*(z4_10_1 - 
	1.0) + 5.0d-1*(x5_10_1 - 1.0)*(x5_10_1 - 1.0) + 5.0d-1*(z5_10_1 - 1.0)*(z5_10_1 
	- 1.0) + 5.0d-1*(x6_10_1 - 1.0)*(x6_10_1 - 1.0) + 5.0d-1*(z6_10_1 - 
	1.0)*(z6_10_1 - 1.0) + 5.0d-1*(x7_10_1 - 1.0)*(x7_10_1 - 1.0) + 5.0d-1*(z7_10_1 
	- 1.0)*(z7_10_1 - 1.0) + 5.0d-1*(x8_10_1 - 1.0)*(x8_10_1 - 1.0) + 
	5.0d-1*(z8_10_1 - 1.0)*(z8_10_1 - 1.0) + 5.0d-1*(x9_10_1 - 1.0)*(x9_10_1 - 1.0) 
	+ 5.0d-1*(z9_10_1 - 1.0)*(z9_10_1 - 1.0) + 5.0d-1*(x1_10_2 - 1.0)*(x1_10_2 - 
	1.0) + 5.0d-1*(z1_10_2 - 1.0)*(z1_10_2 - 1.0) + 5.0d-1*(x2_10_2 - 1.0)*(x2_10_2 
	- 1.0) + 5.0d-1*(z2_10_2 - 1.0)*(z2_10_2 - 1.0) + 5.0d-1*(x3_10_2 - 
	1.0)*(x3_10_2 - 1.0) + 5.0d-1*(z3_10_2 - 1.0)*(z3_10_2 - 1.0) + 5.0d-1*(x4_10_2 
	- 1.0)*(x4_10_2 - 1.0) + 5.0d-1*(z4_10_2 - 1.0)*(z4_10_2 - 1.0) + 
	5.0d-1*(x5_10_2 - 1.0)*(x5_10_2 - 1.0) + 5.0d-1*(z5_10_2 - 1.0)*(z5_10_2 - 1.0) 
	+ 5.0d-1*(x6_10_2 - 1.0)*(x6_10_2 - 1.0) + 5.0d-1*(z6_10_2 - 1.0)*(z6_10_2 - 
	1.0) + 5.0d-1*(x7_10_2 - 1.0)*(x7_10_2 - 1.0) + 5.0d-1*(z7_10_2 - 1.0)*(z7_10_2 
	- 1.0) + 5.0d-1*(x8_10_2 - 1.0)*(x8_10_2 - 1.0) + 5.0d-1*(z8_10_2 - 
	1.0)*(z8_10_2 - 1.0) + 5.0d-1*(x9_10_2 - 1.0)*(x9_10_2 - 1.0) + 5.0d-1*(z9_10_2 
	- 1.0)*(z9_10_2 - 1.0) + 5.0d-1*(x1_10_3 - 1.0)*(x1_10_3 - 1.0) + 
	5.0d-1*(z1_10_3 - 1.0)*(z1_10_3 - 1.0) + 5.0d-1*(x2_10_3 - 1.0)*(x2_10_3 - 1.0) 
	+ 5.0d-1*(z2_10_3 - 1.0)*(z2_10_3 - 1.0) + 5.0d-1*(x3_10_3 - 1.0)*(x3_10_3 - 
	1.0) + 5.0d-1*(z3_10_3 - 1.0)*(z3_10_3 - 1.0) + 5.0d-1*(x4_10_3 - 1.0)*(x4_10_3 
	- 1.0) + 5.0d-1*(z4_10_3 - 1.0)*(z4_10_3 - 1.0) + 5.0d-1*(x5_10_3 - 
	1.0)*(x5_10_3 - 1.0) + 5.0d-1*(z5_10_3 - 1.0)*(z5_10_3 - 1.0) + 5.0d-1*(x6_10_3 
	- 1.0)*(x6_10_3 - 1.0) + 5.0d-1*(z6_10_3 - 1.0)*(z6_10_3 - 1.0) + 
	5.0d-1*(x7_10_3 - 1.0)*(x7_10_3 - 1.0) + 5.0d-1*(z7_10_3 - 1.0)*(z7_10_3 - 1.0) 
	+ 5.0d-1*(x8_10_3 - 1.0)*(x8_10_3 - 1.0) + 5.0d-1*(z8_10_3 - 1.0)*(z8_10_3 - 
	1.0) + 5.0d-1*(x9_10_3 - 1.0)*(x9_10_3 - 1.0) + 5.0d-1*(z9_10_3 - 1.0)*(z9_10_3 
	- 1.0) + 5.0d-1*(x1_10_4 - 1.0)*(x1_10_4 - 1.0) + 5.0d-1*(z1_10_4 - 
	1.0)*(z1_10_4 - 1.0) + 5.0d-1*(x2_10_4 - 1.0)*(x2_10_4 - 1.0) + 5.0d-1*(z2_10_4 
	- 1.0)*(z2_10_4 - 1.0) + 5.0d-1*(x3_10_4 - 1.0)*(x3_10_4 - 1.0) + 
	5.0d-1*(z3_10_4 - 1.0)*(z3_10_4 - 1.0) + 5.0d-1*(x4_10_4 - 1.0)*(x4_10_4 - 1.0) 
	+ 5.0d-1*(z4_10_4 - 1.0)*(z4_10_4 - 1.0) + 5.0d-1*(x5_10_4 - 1.0)*(x5_10_4 - 
	1.0) + 5.0d-1*(z5_10_4 - 1.0)*(z5_10_4 - 1.0) + 5.0d-1*(x6_10_4 - 1.0)*(x6_10_4 
	- 1.0) + 5.0d-1*(z6_10_4 - 1.0)*(z6_10_4 - 1.0) + 5.0d-1*(x7_10_4 - 
	1.0)*(x7_10_4 - 1.0) + 5.0d-1*(z7_10_4 - 1.0)*(z7_10_4 - 1.0) + 5.0d-1*(x8_10_4 
	- 1.0)*(x8_10_4 - 1.0) + 5.0d-1*(z8_10_4 - 1.0)*(z8_10_4 - 1.0) + 
	5.0d-1*(x9_10_4 - 1.0)*(x9_10_4 - 1.0) + 5.0d-1*(z9_10_4 - 1.0)*(z9_10_4 - 1.0) 
	+ 5.0d-1*(x1_10_5 - 1.0)*(x1_10_5 - 1.0) + 5.0d-1*(z1_10_5 - 1.0)*(z1_10_5 - 
	1.0) + 5.0d-1*(x2_10_5 - 1.0)*(x2_10_5 - 1.0) + 5.0d-1*(z2_10_5 - 1.0)*(z2_10_5 
	- 1.0) + 5.0d-1*(x3_10_5 - 1.0)*(x3_10_5 - 1.0) + 5.0d-1*(z3_10_5 - 
	1.0)*(z3_10_5 - 1.0) + 5.0d-1*(x4_10_5 - 1.0)*(x4_10_5 - 1.0) + 5.0d-1*(z4_10_5 
	- 1.0)*(z4_10_5 - 1.0) + 5.0d-1*(x5_10_5 - 1.0)*(x5_10_5 - 1.0) + 
	5.0d-1*(z5_10_5 - 1.0)*(z5_10_5 - 1.0) + 5.0d-1*(x6_10_5 - 1.0)*(x6_10_5 - 1.0) 
	+ 5.0d-1*(z6_10_5 - 1.0)*(z6_10_5 - 1.0) + 5.0d-1*(x7_10_5 - 1.0)*(x7_10_5 - 
	1.0) + 5.0d-1*(z7_10_5 - 1.0)*(z7_10_5 - 1.0) + 5.0d-1*(x8_10_5 - 1.0)*(x8_10_5 
	- 1.0) + 5.0d-1*(z8_10_5 - 1.0)*(z8_10_5 - 1.0) + 5.0d-1*(x9_10_5 - 
	1.0)*(x9_10_5 - 1.0) + 5.0d-1*(z9_10_5 - 1.0)*(z9_10_5 - 1.0) + 5.0d-1*(x1_10_6 
	- 1.0)*(x1_10_6 - 1.0) + 5.0d-1*(z1_10_6 - 1.0)*(z1_10_6 - 1.0) + 
	5.0d-1*(x2_10_6 - 1.0)*(x2_10_6 - 1.0) + 5.0d-1*(z2_10_6 - 1.0)*(z2_10_6 - 1.0) 
	+ 5.0d-1*(x3_10_6 - 1.0)*(x3_10_6 - 1.0) + 5.0d-1*(z3_10_6 - 1.0)*(z3_10_6 - 
	1.0) + 5.0d-1*(x4_10_6 - 1.0)*(x4_10_6 - 1.0) + 5.0d-1*(z4_10_6 - 1.0)*(z4_10_6 
	- 1.0) + 5.0d-1*(x5_10_6 - 1.0)*(x5_10_6 - 1.0) + 5.0d-1*(z5_10_6 - 
	1.0)*(z5_10_6 - 1.0) + 5.0d-1*(x6_10_6 - 1.0)*(x6_10_6 - 1.0) + 5.0d-1*(z6_10_6 
	- 1.0)*(z6_10_6 - 1.0) + 5.0d-1*(x7_10_6 - 1.0)*(x7_10_6 - 1.0) + 
	5.0d-1*(z7_10_6 - 1.0)*(z7_10_6 - 1.0) + 5.0d-1*(x8_10_6 - 1.0)*(x8_10_6 - 1.0) 
	+ 5.0d-1*(z8_10_6 - 1.0)*(z8_10_6 - 1.0) + 5.0d-1*(x9_10_6 - 1.0)*(x9_10_6 - 
	1.0) + 5.0d-1*(z9_10_6 - 1.0)*(z9_10_6 - 1.0) + 5.0d-1*(x1_10_7 - 1.0)*(x1_10_7 
	- 1.0) + 5.0d-1*(z1_10_7 - 1.0)*(z1_10_7 - 1.0) + 5.0d-1*(x2_10_7 - 
	1.0)*(x2_10_7 - 1.0) + 5.0d-1*(z2_10_7 - 1.0)*(z2_10_7 - 1.0) + 5.0d-1*(x3_10_7 
	- 1.0)*(x3_10_7 - 1.0) + 5.0d-1*(z3_10_7 - 1.0)*(z3_10_7 - 1.0) + 
	5.0d-1*(x4_10_7 - 1.0)*(x4_10_7 - 1.0) + 5.0d-1*(z4_10_7 - 1.0)*(z4_10_7 - 1.0) 
	+ 5.0d-1*(x5_10_7 - 1.0)*(x5_10_7 - 1.0) + 5.0d-1*(z5_10_7 - 1.0)*(z5_10_7 - 
	1.0) + 5.0d-1*(x6_10_7 - 1.0)*(x6_10_7 - 1.0) + 5.0d-1*(z6_10_7 - 1.0)*(z6_10_7 
	- 1.0) + 5.0d-1*(x7_10_7 - 1.0)*(x7_10_7 - 1.0) + 5.0d-1*(z7_10_7 - 
	1.0)*(z7_10_7 - 1.0) + 5.0d-1*(x8_10_7 - 1.0)*(x8_10_7 - 1.0) + 5.0d-1*(z8_10_7 
	- 1.0)*(z8_10_7 - 1.0) + 5.0d-1*(x9_10_7 - 1.0)*(x9_10_7 - 1.0) + 
	5.0d-1*(z9_10_7 - 1.0)*(z9_10_7 - 1.0) + 5.0d-1*(x1_10_8 - 1.0)*(x1_10_8 - 1.0) 
	+ 5.0d-1*(z1_10_8 - 1.0)*(z1_10_8 - 1.0) + 5.0d-1*(x2_10_8 - 1.0)*(x2_10_8 - 
	1.0) + 5.0d-1*(z2_10_8 - 1.0)*(z2_10_8 - 1.0) + 5.0d-1*(x3_10_8 - 1.0)*(x3_10_8 
	- 1.0) + 5.0d-1*(z3_10_8 - 1.0)*(z3_10_8 - 1.0) + 5.0d-1*(x4_10_8 - 
	1.0)*(x4_10_8 - 1.0) + 5.0d-1*(z4_10_8 - 1.0)*(z4_10_8 - 1.0) + 5.0d-1*(x5_10_8 
	- 1.0)*(x5_10_8 - 1.0) + 5.0d-1*(z5_10_8 - 1.0)*(z5_10_8 - 1.0) + 
	5.0d-1*(x6_10_8 - 1.0)*(x6_10_8 - 1.0) + 5.0d-1*(z6_10_8 - 1.0)*(z6_10_8 - 1.0) 
	+ 5.0d-1*(x7_10_8 - 1.0)*(x7_10_8 - 1.0) + 5.0d-1*(z7_10_8 - 1.0)*(z7_10_8 - 
	1.0) + 5.0d-1*(x8_10_8 - 1.0)*(x8_10_8 - 1.0) + 5.0d-1*(z8_10_8 - 1.0)*(z8_10_8 
	- 1.0) + 5.0d-1*(x9_10_8 - 1.0)*(x9_10_8 - 1.0) + 5.0d-1*(z9_10_8 - 
	1.0)*(z9_10_8 - 1.0) + 5.0d-1*(x1_10_9 - 1.0)*(x1_10_9 - 1.0) + 5.0d-1*(z1_10_9 
	- 1.0)*(z1_10_9 - 1.0) + 5.0d-1*(x2_10_9 - 1.0)*(x2_10_9 - 1.0) + 
	5.0d-1*(z2_10_9 - 1.0)*(z2_10_9 - 1.0) + 5.0d-1*(x3_10_9 - 1.0)*(x3_10_9 - 1.0) 
	+ 5.0d-1*(z3_10_9 - 1.0)*(z3_10_9 - 1.0) + 5.0d-1*(x4_10_9 - 1.0)*(x4_10_9 - 
	1.0) + 5.0d-1*(z4_10_9 - 1.0)*(z4_10_9 - 1.0) + 5.0d-1*(x5_10_9 - 1.0)*(x5_10_9 
	- 1.0) + 5.0d-1*(z5_10_9 - 1.0)*(z5_10_9 - 1.0) + 5.0d-1*(x6_10_9 - 
	1.0)*(x6_10_9 - 1.0) + 5.0d-1*(z6_10_9 - 1.0)*(z6_10_9 - 1.0) + 5.0d-1*(x7_10_9 
	- 1.0)*(x7_10_9 - 1.0) + 5.0d-1*(z7_10_9 - 1.0)*(z7_10_9 - 1.0) + 
	5.0d-1*(x8_10_9 - 1.0)*(x8_10_9 - 1.0) + 5.0d-1*(z8_10_9 - 1.0)*(z8_10_9 - 1.0) 
	+ 5.0d-1*(x9_10_9 - 1.0)*(x9_10_9 - 1.0) + 5.0d-1*(z9_10_9 - 1.0)*(z9_10_9 - 
	1.0) + 5.0d-1*(x1_1_10 - 1.0)*(x1_1_10 - 1.0) + 5.0d-1*(y1_1_10 - 1.0)*(y1_1_10 
	- 1.0) + 5.0d-1*(x2_1_10 - 1.0)*(x2_1_10 - 1.0) + 5.0d-1*(y2_1_10 - 
	1.0)*(y2_1_10 - 1.0) + 5.0d-1*(x3_1_10 - 1.0)*(x3_1_10 - 1.0) + 5.0d-1*(y3_1_10 
	- 1.0)*(y3_1_10 - 1.0) + 5.0d-1*(x4_1_10 - 1.0)*(x4_1_10 - 1.0) + 
	5.0d-1*(y4_1_10 - 1.0)*(y4_1_10 - 1.0) + 5.0d-1*(x5_1_10 - 1.0)*(x5_1_10 - 1.0) 
	+ 5.0d-1*(y5_1_10 - 1.0)*(y5_1_10 - 1.0) + 5.0d-1*(x6_1_10 - 1.0)*(x6_1_10 - 
	1.0) + 5.0d-1*(y6_1_10 - 1.0)*(y6_1_10 - 1.0) + 5.0d-1*(x7_1_10 - 1.0)*(x7_1_10 
	- 1.0) + 5.0d-1*(y7_1_10 - 1.0)*(y7_1_10 - 1.0) + 5.0d-1*(x8_1_10 - 
	1.0)*(x8_1_10 - 1.0) + 5.0d-1*(y8_1_10 - 1.0)*(y8_1_10 - 1.0) + 5.0d-1*(x9_1_10 
	- 1.0)*(x9_1_10 - 1.0) + 5.0d-1*(y9_1_10 - 1.0)*(y9_1_10 - 1.0) + 
	5.0d-1*(x1_2_10 - 1.0)*(x1_2_10 - 1.0) + 5.0d-1*(y1_2_10 - 1.0)*(y1_2_10 - 1.0) 
	+ 5.0d-1*(x2_2_10 - 1.0)*(x2_2_10 - 1.0) + 5.0d-1*(y2_2_10 - 1.0)*(y2_2_10 - 
	1.0) + 5.0d-1*(x3_2_10 - 1.0)*(x3_2_10 - 1.0) + 5.0d-1*(y3_2_10 - 1.0)*(y3_2_10 
	- 1.0) + 5.0d-1*(x4_2_10 - 1.0)*(x4_2_10 - 1.0) + 5.0d-1*(y4_2_10 - 
	1.0)*(y4_2_10 - 1.0) + 5.0d-1*(x5_2_10 - 1.0)*(x5_2_10 - 1.0) + 5.0d-1*(y5_2_10 
	- 1.0)*(y5_2_10 - 1.0) + 5.0d-1*(x6_2_10 - 1.0)*(x6_2_10 - 1.0) + 
	5.0d-1*(y6_2_10 - 1.0)*(y6_2_10 - 1.0) + 5.0d-1*(x7_2_10 - 1.0)*(x7_2_10 - 1.0) 
	+ 5.0d-1*(y7_2_10 - 1.0)*(y7_2_10 - 1.0) + 5.0d-1*(x8_2_10 - 1.0)*(x8_2_10 - 
	1.0) + 5.0d-1*(y8_2_10 - 1.0)*(y8_2_10 - 1.0) + 5.0d-1*(x9_2_10 - 1.0)*(x9_2_10 
	- 1.0) + 5.0d-1*(y9_2_10 - 1.0)*(y9_2_10 - 1.0) + 5.0d-1*(x1_3_10 - 
	1.0)*(x1_3_10 - 1.0) + 5.0d-1*(y1_3_10 - 1.0)*(y1_3_10 - 1.0) + 5.0d-1*(x2_3_10 
	- 1.0)*(x2_3_10 - 1.0) + 5.0d-1*(y2_3_10 - 1.0)*(y2_3_10 - 1.0) + 
	5.0d-1*(x3_3_10 - 1.0)*(x3_3_10 - 1.0) + 5.0d-1*(y3_3_10 - 1.0)*(y3_3_10 - 1.0) 
	+ 5.0d-1*(x4_3_10 - 1.0)*(x4_3_10 - 1.0) + 5.0d-1*(y4_3_10 - 1.0)*(y4_3_10 - 
	1.0) + 5.0d-1*(x5_3_10 - 1.0)*(x5_3_10 - 1.0) + 5.0d-1*(y5_3_10 - 1.0)*(y5_3_10 
	- 1.0) + 5.0d-1*(x6_3_10 - 1.0)*(x6_3_10 - 1.0) + 5.0d-1*(y6_3_10 - 
	1.0)*(y6_3_10 - 1.0) + 5.0d-1*(x7_3_10 - 1.0)*(x7_3_10 - 1.0) + 5.0d-1*(y7_3_10 
	- 1.0)*(y7_3_10 - 1.0) + 5.0d-1*(x8_3_10 - 1.0)*(x8_3_10 - 1.0) + 
	5.0d-1*(y8_3_10 - 1.0)*(y8_3_10 - 1.0) + 5.0d-1*(x9_3_10 - 1.0)*(x9_3_10 - 1.0) 
	+ 5.0d-1*(y9_3_10 - 1.0)*(y9_3_10 - 1.0) + 5.0d-1*(x1_4_10 - 1.0)*(x1_4_10 - 
	1.0) + 5.0d-1*(y1_4_10 - 1.0)*(y1_4_10 - 1.0) + 5.0d-1*(x2_4_10 - 1.0)*(x2_4_10 
	- 1.0) + 5.0d-1*(y2_4_10 - 1.0)*(y2_4_10 - 1.0) + 5.0d-1*(x3_4_10 - 
	1.0)*(x3_4_10 - 1.0) + 5.0d-1*(y3_4_10 - 1.0)*(y3_4_10 - 1.0) + 5.0d-1*(x4_4_10 
	- 1.0)*(x4_4_10 - 1.0) + 5.0d-1*(y4_4_10 - 1.0)*(y4_4_10 - 1.0) + 
	5.0d-1*(x5_4_10 - 1.0)*(x5_4_10 - 1.0) + 5.0d-1*(y5_4_10 - 1.0)*(y5_4_10 - 1.0) 
	+ 5.0d-1*(x6_4_10 - 1.0)*(x6_4_10 - 1.0) + 5.0d-1*(y6_4_10 - 1.0)*(y6_4_10 - 
	1.0) + 5.0d-1*(x7_4_10 - 1.0)*(x7_4_10 - 1.0) + 5.0d-1*(y7_4_10 - 1.0)*(y7_4_10 
	- 1.0) + 5.0d-1*(x8_4_10 - 1.0)*(x8_4_10 - 1.0) + 5.0d-1*(y8_4_10 - 
	1.0)*(y8_4_10 - 1.0) + 5.0d-1*(x9_4_10 - 1.0)*(x9_4_10 - 1.0) + 5.0d-1*(y9_4_10 
	- 1.0)*(y9_4_10 - 1.0) + 5.0d-1*(x1_5_10 - 1.0)*(x1_5_10 - 1.0) + 
	5.0d-1*(y1_5_10 - 1.0)*(y1_5_10 - 1.0) + 5.0d-1*(x2_5_10 - 1.0)*(x2_5_10 - 1.0) 
	+ 5.0d-1*(y2_5_10 - 1.0)*(y2_5_10 - 1.0) + 5.0d-1*(x3_5_10 - 1.0)*(x3_5_10 - 
	1.0) + 5.0d-1*(y3_5_10 - 1.0)*(y3_5_10 - 1.0) + 5.0d-1*(x4_5_10 - 1.0)*(x4_5_10 
	- 1.0) + 5.0d-1*(y4_5_10 - 1.0)*(y4_5_10 - 1.0) + 5.0d-1*(x5_5_10 - 
	1.0)*(x5_5_10 - 1.0) + 5.0d-1*(y5_5_10 - 1.0)*(y5_5_10 - 1.0) + 5.0d-1*(x6_5_10 
	- 1.0)*(x6_5_10 - 1.0) + 5.0d-1*(y6_5_10 - 1.0)*(y6_5_10 - 1.0) + 
	5.0d-1*(x7_5_10 - 1.0)*(x7_5_10 - 1.0) + 5.0d-1*(y7_5_10 - 1.0)*(y7_5_10 - 1.0) 
	+ 5.0d-1*(x8_5_10 - 1.0)*(x8_5_10 - 1.0) + 5.0d-1*(y8_5_10 - 1.0)*(y8_5_10 - 
	1.0) + 5.0d-1*(x9_5_10 - 1.0)*(x9_5_10 - 1.0) + 5.0d-1*(y9_5_10 - 1.0)*(y9_5_10 
	- 1.0) + 5.0d-1*(x1_6_10 - 1.0)*(x1_6_10 - 1.0) + 5.0d-1*(y1_6_10 - 
	1.0)*(y1_6_10 - 1.0) + 5.0d-1*(x2_6_10 - 1.0)*(x2_6_10 - 1.0) + 5.0d-1*(y2_6_10 
	- 1.0)*(y2_6_10 - 1.0) + 5.0d-1*(x3_6_10 - 1.0)*(x3_6_10 - 1.0) + 
	5.0d-1*(y3_6_10 - 1.0)*(y3_6_10 - 1.0) + 5.0d-1*(x4_6_10 - 1.0)*(x4_6_10 - 1.0) 
	+ 5.0d-1*(y4_6_10 - 1.0)*(y4_6_10 - 1.0) + 5.0d-1*(x5_6_10 - 1.0)*(x5_6_10 - 
	1.0) + 5.0d-1*(y5_6_10 - 1.0)*(y5_6_10 - 1.0) + 5.0d-1*(x6_6_10 - 1.0)*(x6_6_10 
	- 1.0) + 5.0d-1*(y6_6_10 - 1.0)*(y6_6_10 - 1.0) + 5.0d-1*(x7_6_10 - 
	1.0)*(x7_6_10 - 1.0) + 5.0d-1*(y7_6_10 - 1.0)*(y7_6_10 - 1.0) + 5.0d-1*(x8_6_10 
	- 1.0)*(x8_6_10 - 1.0) + 5.0d-1*(y8_6_10 - 1.0)*(y8_6_10 - 1.0) + 
	5.0d-1*(x9_6_10 - 1.0)*(x9_6_10 - 1.0) + 5.0d-1*(y9_6_10 - 1.0)*(y9_6_10 - 1.0) 
	+ 5.0d-1*(x1_7_10 - 1.0)*(x1_7_10 - 1.0) + 5.0d-1*(y1_7_10 - 1.0)*(y1_7_10 - 
	1.0) + 5.0d-1*(x2_7_10 - 1.0)*(x2_7_10 - 1.0) + 5.0d-1*(y2_7_10 - 1.0)*(y2_7_10 
	- 1.0) + 5.0d-1*(x3_7_10 - 1.0)*(x3_7_10 - 1.0) + 5.0d-1*(y3_7_10 - 
	1.0)*(y3_7_10 - 1.0) + 5.0d-1*(x4_7_10 - 1.0)*(x4_7_10 - 1.0) + 5.0d-1*(y4_7_10 
	- 1.0)*(y4_7_10 - 1.0) + 5.0d-1*(x5_7_10 - 1.0)*(x5_7_10 - 1.0) + 
	5.0d-1*(y5_7_10 - 1.0)*(y5_7_10 - 1.0) + 5.0d-1*(x6_7_10 - 1.0)*(x6_7_10 - 1.0) 
	+ 5.0d-1*(y6_7_10 - 1.0)*(y6_7_10 - 1.0) + 5.0d-1*(x7_7_10 - 1.0)*(x7_7_10 - 
	1.0) + 5.0d-1*(y7_7_10 - 1.0)*(y7_7_10 - 1.0) + 5.0d-1*(x8_7_10 - 1.0)*(x8_7_10 
	- 1.0) + 5.0d-1*(y8_7_10 - 1.0)*(y8_7_10 - 1.0) + 5.0d-1*(x9_7_10 - 
	1.0)*(x9_7_10 - 1.0) + 5.0d-1*(y9_7_10 - 1.0)*(y9_7_10 - 1.0) + 5.0d-1*(x1_8_10 
	- 1.0)*(x1_8_10 - 1.0) + 5.0d-1*(y1_8_10 - 1.0)*(y1_8_10 - 1.0) + 
	5.0d-1*(x2_8_10 - 1.0)*(x2_8_10 - 1.0) + 5.0d-1*(y2_8_10 - 1.0)*(y2_8_10 - 1.0) 
	+ 5.0d-1*(x3_8_10 - 1.0)*(x3_8_10 - 1.0) + 5.0d-1*(y3_8_10 - 1.0)*(y3_8_10 - 
	1.0) + 5.0d-1*(x4_8_10 - 1.0)*(x4_8_10 - 1.0) + 5.0d-1*(y4_8_10 - 1.0)*(y4_8_10 
	- 1.0) + 5.0d-1*(x5_8_10 - 1.0)*(x5_8_10 - 1.0) + 5.0d-1*(y5_8_10 - 
	1.0)*(y5_8_10 - 1.0) + 5.0d-1*(x6_8_10 - 1.0)*(x6_8_10 - 1.0) + 5.0d-1*(y6_8_10 
	- 1.0)*(y6_8_10 - 1.0) + 5.0d-1*(x7_8_10 - 1.0)*(x7_8_10 - 1.0) + 
	5.0d-1*(y7_8_10 - 1.0)*(y7_8_10 - 1.0) + 5.0d-1*(x8_8_10 - 1.0)*(x8_8_10 - 1.0) 
	+ 5.0d-1*(y8_8_10 - 1.0)*(y8_8_10 - 1.0) + 5.0d-1*(x9_8_10 - 1.0)*(x9_8_10 - 
	1.0) + 5.0d-1*(y9_8_10 - 1.0)*(y9_8_10 - 1.0) + 5.0d-1*(x1_9_10 - 1.0)*(x1_9_10 
	- 1.0) + 5.0d-1*(y1_9_10 - 1.0)*(y1_9_10 - 1.0) + 5.0d-1*(x2_9_10 - 
	1.0)*(x2_9_10 - 1.0) + 5.0d-1*(y2_9_10 - 1.0)*(y2_9_10 - 1.0) + 5.0d-1*(x3_9_10 
	- 1.0)*(x3_9_10 - 1.0) + 5.0d-1*(y3_9_10 - 1.0)*(y3_9_10 - 1.0) + 
	5.0d-1*(x4_9_10 - 1.0)*(x4_9_10 - 1.0) + 5.0d-1*(y4_9_10 - 1.0)*(y4_9_10 - 1.0) 
	+ 5.0d-1*(x5_9_10 - 1.0)*(x5_9_10 - 1.0) + 5.0d-1*(y5_9_10 - 1.0)*(y5_9_10 - 
	1.0) + 5.0d-1*(x6_9_10 - 1.0)*(x6_9_10 - 1.0) + 5.0d-1*(y6_9_10 - 1.0)*(y6_9_10 
	- 1.0) + 5.0d-1*(x7_9_10 - 1.0)*(x7_9_10 - 1.0) + 5.0d-1*(y7_9_10 - 
	1.0)*(y7_9_10 - 1.0) + 5.0d-1*(x8_9_10 - 1.0)*(x8_9_10 - 1.0) + 5.0d-1*(y8_9_10 
	- 1.0)*(y8_9_10 - 1.0) + 5.0d-1*(x9_9_10 - 1.0)*(x9_9_10 - 1.0) + 
	5.0d-1*(y9_9_10 - 1.0)*(y9_9_10 - 1.0);

subject to v1_1_1:
	x1_1_1 + y1_1_1 + z1_1_1 + y0_1_1 + z0_1_1 + x1_0_1 + z1_0_1 + x1_1_0 + y1_1_0 
	- 1.0 = 0;
subject to v2_1_1:
	-x1_1_1 + x2_1_1 + y2_1_1 + z2_1_1 + x2_0_1 + z2_0_1 + x2_1_0 + y2_1_0 - 1.0 = 
	0;
subject to v3_1_1:
	-x2_1_1 + x3_1_1 + y3_1_1 + z3_1_1 + x3_0_1 + z3_0_1 + x3_1_0 + y3_1_0 - 1.0 = 
	0;
subject to v4_1_1:
	-x3_1_1 + x4_1_1 + y4_1_1 + z4_1_1 + x4_0_1 + z4_0_1 + x4_1_0 + y4_1_0 - 1.0 = 
	0;
subject to v5_1_1:
	-x4_1_1 + x5_1_1 + y5_1_1 + z5_1_1 + x5_0_1 + z5_0_1 + x5_1_0 + y5_1_0 - 1.0 = 
	0;
subject to v6_1_1:
	-x5_1_1 + x6_1_1 + y6_1_1 + z6_1_1 + x6_0_1 + z6_0_1 + x6_1_0 + y6_1_0 - 1.0 = 
	0;
subject to v7_1_1:
	-x6_1_1 + x7_1_1 + y7_1_1 + z7_1_1 + x7_0_1 + z7_0_1 + x7_1_0 + y7_1_0 - 1.0 = 
	0;
subject to v8_1_1:
	-x7_1_1 + x8_1_1 + y8_1_1 + z8_1_1 + x8_0_1 + z8_0_1 + x8_1_0 + y8_1_0 - 1.0 = 
	0;
subject to v9_1_1:
	-x8_1_1 + x9_1_1 + y9_1_1 + z9_1_1 + x9_0_1 + z9_0_1 + x9_1_0 + y9_1_0 - 1.0 = 
	0;
subject to v10_1_1:
	-x9_1_1 + y10_1_1 + z10_1_1 + y11_1_1 + z11_1_1 + x10_0_1 + z10_0_1 + x10_1_0 + 
	y10_1_0 - 1.0 = 0;
subject to v1_2_1:
	-y1_1_1 + x1_2_1 + y1_2_1 + z1_2_1 + y0_2_1 + z0_2_1 + x1_2_0 + y1_2_0 - 1.0 = 
	0;
subject to v2_2_1:
	-y2_1_1 - x1_2_1 + x2_2_1 + y2_2_1 + z2_2_1 + x2_2_0 + y2_2_0 - 1.0 = 0;
subject to v3_2_1:
	-y3_1_1 - x2_2_1 + x3_2_1 + y3_2_1 + z3_2_1 + x3_2_0 + y3_2_0 - 1.0 = 0;
subject to v4_2_1:
	-y4_1_1 - x3_2_1 + x4_2_1 + y4_2_1 + z4_2_1 + x4_2_0 + y4_2_0 - 1.0 = 0;
subject to v5_2_1:
	-y5_1_1 - x4_2_1 + x5_2_1 + y5_2_1 + z5_2_1 + x5_2_0 + y5_2_0 - 1.0 = 0;
subject to v6_2_1:
	-y6_1_1 - x5_2_1 + x6_2_1 + y6_2_1 + z6_2_1 + x6_2_0 + y6_2_0 - 1.0 = 0;
subject to v7_2_1:
	-y7_1_1 - x6_2_1 + x7_2_1 + y7_2_1 + z7_2_1 + x7_2_0 + y7_2_0 - 1.0 = 0;
subject to v8_2_1:
	-y8_1_1 - x7_2_1 + x8_2_1 + y8_2_1 + z8_2_1 + x8_2_0 + y8_2_0 - 1.0 = 0;
subject to v9_2_1:
	-y9_1_1 - x8_2_1 + x9_2_1 + y9_2_1 + z9_2_1 + x9_2_0 + y9_2_0 - 1.0 = 0;
subject to v10_2_1:
	-x9_2_1 - y10_1_1 + y10_2_1 + z10_2_1 + y11_2_1 + z11_2_1 + x10_2_0 + y10_2_0 - 
	1.0 = 0;
subject to v1_3_1:
	-y1_2_1 + x1_3_1 + y1_3_1 + z1_3_1 + y0_3_1 + z0_3_1 + x1_3_0 + y1_3_0 - 1.0 = 
	0;
subject to v2_3_1:
	-y2_2_1 - x1_3_1 + x2_3_1 + y2_3_1 + z2_3_1 + x2_3_0 + y2_3_0 - 1.0 = 0;
subject to v3_3_1:
	-y3_2_1 - x2_3_1 + x3_3_1 + y3_3_1 + z3_3_1 + x3_3_0 + y3_3_0 - 1.0 = 0;
subject to v4_3_1:
	-y4_2_1 - x3_3_1 + x4_3_1 + y4_3_1 + z4_3_1 + x4_3_0 + y4_3_0 - 1.0 = 0;
subject to v5_3_1:
	-y5_2_1 - x4_3_1 + x5_3_1 + y5_3_1 + z5_3_1 + x5_3_0 + y5_3_0 - 1.0 = 0;
subject to v6_3_1:
	-y6_2_1 - x5_3_1 + x6_3_1 + y6_3_1 + z6_3_1 + x6_3_0 + y6_3_0 - 1.0 = 0;
subject to v7_3_1:
	-y7_2_1 - x6_3_1 + x7_3_1 + y7_3_1 + z7_3_1 + x7_3_0 + y7_3_0 - 1.0 = 0;
subject to v8_3_1:
	-y8_2_1 - x7_3_1 + x8_3_1 + y8_3_1 + z8_3_1 + x8_3_0 + y8_3_0 - 1.0 = 0;
subject to v9_3_1:
	-y9_2_1 - x8_3_1 + x9_3_1 + y9_3_1 + z9_3_1 + x9_3_0 + y9_3_0 - 1.0 = 0;
subject to v10_3_1:
	-x9_3_1 - y10_2_1 + y10_3_1 + z10_3_1 + y11_3_1 + z11_3_1 + x10_3_0 + y10_3_0 - 
	1.0 = 0;
subject to v1_4_1:
	-y1_3_1 + x1_4_1 + y1_4_1 + z1_4_1 + y0_4_1 + z0_4_1 + x1_4_0 + y1_4_0 - 1.0 = 
	0;
subject to v2_4_1:
	-y2_3_1 - x1_4_1 + x2_4_1 + y2_4_1 + z2_4_1 + x2_4_0 + y2_4_0 - 1.0 = 0;
subject to v3_4_1:
	-y3_3_1 - x2_4_1 + x3_4_1 + y3_4_1 + z3_4_1 + x3_4_0 + y3_4_0 - 1.0 = 0;
subject to v4_4_1:
	-y4_3_1 - x3_4_1 + x4_4_1 + y4_4_1 + z4_4_1 + x4_4_0 + y4_4_0 - 1.0 = 0;
subject to v5_4_1:
	-y5_3_1 - x4_4_1 + x5_4_1 + y5_4_1 + z5_4_1 + x5_4_0 + y5_4_0 - 1.0 = 0;
subject to v6_4_1:
	-y6_3_1 - x5_4_1 + x6_4_1 + y6_4_1 + z6_4_1 + x6_4_0 + y6_4_0 - 1.0 = 0;
subject to v7_4_1:
	-y7_3_1 - x6_4_1 + x7_4_1 + y7_4_1 + z7_4_1 + x7_4_0 + y7_4_0 - 1.0 = 0;
subject to v8_4_1:
	-y8_3_1 - x7_4_1 + x8_4_1 + y8_4_1 + z8_4_1 + x8_4_0 + y8_4_0 - 1.0 = 0;
subject to v9_4_1:
	-y9_3_1 - x8_4_1 + x9_4_1 + y9_4_1 + z9_4_1 + x9_4_0 + y9_4_0 - 1.0 = 0;
subject to v10_4_1:
	-x9_4_1 - y10_3_1 + y10_4_1 + z10_4_1 + y11_4_1 + z11_4_1 + x10_4_0 + y10_4_0 - 
	1.0 = 0;
subject to v1_5_1:
	-y1_4_1 + x1_5_1 + y1_5_1 + z1_5_1 + y0_5_1 + z0_5_1 + x1_5_0 + y1_5_0 - 1.0 = 
	0;
subject to v2_5_1:
	-y2_4_1 - x1_5_1 + x2_5_1 + y2_5_1 + z2_5_1 + x2_5_0 + y2_5_0 - 1.0 = 0;
subject to v3_5_1:
	-y3_4_1 - x2_5_1 + x3_5_1 + y3_5_1 + z3_5_1 + x3_5_0 + y3_5_0 - 1.0 = 0;
subject to v4_5_1:
	-y4_4_1 - x3_5_1 + x4_5_1 + y4_5_1 + z4_5_1 + x4_5_0 + y4_5_0 - 1.0 = 0;
subject to v5_5_1:
	-y5_4_1 - x4_5_1 + x5_5_1 + y5_5_1 + z5_5_1 + x5_5_0 + y5_5_0 - 1.0 = 0;
subject to v6_5_1:
	-y6_4_1 - x5_5_1 + x6_5_1 + y6_5_1 + z6_5_1 + x6_5_0 + y6_5_0 - 1.0 = 0;
subject to v7_5_1:
	-y7_4_1 - x6_5_1 + x7_5_1 + y7_5_1 + z7_5_1 + x7_5_0 + y7_5_0 - 1.0 = 0;
subject to v8_5_1:
	-y8_4_1 - x7_5_1 + x8_5_1 + y8_5_1 + z8_5_1 + x8_5_0 + y8_5_0 - 1.0 = 0;
subject to v9_5_1:
	-y9_4_1 - x8_5_1 + x9_5_1 + y9_5_1 + z9_5_1 + x9_5_0 + y9_5_0 - 1.0 = 0;
subject to v10_5_1:
	-x9_5_1 - y10_4_1 + y10_5_1 + z10_5_1 + y11_5_1 + z11_5_1 + x10_5_0 + y10_5_0 - 
	1.0 = 0;
subject to v1_6_1:
	-y1_5_1 + x1_6_1 + y1_6_1 + z1_6_1 + y0_6_1 + z0_6_1 + x1_6_0 + y1_6_0 - 1.0 = 
	0;
subject to v2_6_1:
	-y2_5_1 - x1_6_1 + x2_6_1 + y2_6_1 + z2_6_1 + x2_6_0 + y2_6_0 - 1.0 = 0;
subject to v3_6_1:
	-y3_5_1 - x2_6_1 + x3_6_1 + y3_6_1 + z3_6_1 + x3_6_0 + y3_6_0 - 1.0 = 0;
subject to v4_6_1:
	-y4_5_1 - x3_6_1 + x4_6_1 + y4_6_1 + z4_6_1 + x4_6_0 + y4_6_0 - 1.0 = 0;
subject to v5_6_1:
	-y5_5_1 - x4_6_1 + x5_6_1 + y5_6_1 + z5_6_1 + x5_6_0 + y5_6_0 - 1.0 = 0;
subject to v6_6_1:
	-y6_5_1 - x5_6_1 + x6_6_1 + y6_6_1 + z6_6_1 + x6_6_0 + y6_6_0 - 1.0 = 0;
subject to v7_6_1:
	-y7_5_1 - x6_6_1 + x7_6_1 + y7_6_1 + z7_6_1 + x7_6_0 + y7_6_0 - 1.0 = 0;
subject to v8_6_1:
	-y8_5_1 - x7_6_1 + x8_6_1 + y8_6_1 + z8_6_1 + x8_6_0 + y8_6_0 - 1.0 = 0;
subject to v9_6_1:
	-y9_5_1 - x8_6_1 + x9_6_1 + y9_6_1 + z9_6_1 + x9_6_0 + y9_6_0 - 1.0 = 0;
subject to v10_6_1:
	-x9_6_1 - y10_5_1 + y10_6_1 + z10_6_1 + y11_6_1 + z11_6_1 + x10_6_0 + y10_6_0 - 
	1.0 = 0;
subject to v1_7_1:
	-y1_6_1 + x1_7_1 + y1_7_1 + z1_7_1 + y0_7_1 + z0_7_1 + x1_7_0 + y1_7_0 - 1.0 = 
	0;
subject to v2_7_1:
	-y2_6_1 - x1_7_1 + x2_7_1 + y2_7_1 + z2_7_1 + x2_7_0 + y2_7_0 - 1.0 = 0;
subject to v3_7_1:
	-y3_6_1 - x2_7_1 + x3_7_1 + y3_7_1 + z3_7_1 + x3_7_0 + y3_7_0 - 1.0 = 0;
subject to v4_7_1:
	-y4_6_1 - x3_7_1 + x4_7_1 + y4_7_1 + z4_7_1 + x4_7_0 + y4_7_0 - 1.0 = 0;
subject to v5_7_1:
	-y5_6_1 - x4_7_1 + x5_7_1 + y5_7_1 + z5_7_1 + x5_7_0 + y5_7_0 - 1.0 = 0;
subject to v6_7_1:
	-y6_6_1 - x5_7_1 + x6_7_1 + y6_7_1 + z6_7_1 + x6_7_0 + y6_7_0 - 1.0 = 0;
subject to v7_7_1:
	-y7_6_1 - x6_7_1 + x7_7_1 + y7_7_1 + z7_7_1 + x7_7_0 + y7_7_0 - 1.0 = 0;
subject to v8_7_1:
	-y8_6_1 - x7_7_1 + x8_7_1 + y8_7_1 + z8_7_1 + x8_7_0 + y8_7_0 - 1.0 = 0;
subject to v9_7_1:
	-y9_6_1 - x8_7_1 + x9_7_1 + y9_7_1 + z9_7_1 + x9_7_0 + y9_7_0 - 1.0 = 0;
subject to v10_7_1:
	-x9_7_1 - y10_6_1 + y10_7_1 + z10_7_1 + y11_7_1 + z11_7_1 + x10_7_0 + y10_7_0 - 
	1.0 = 0;
subject to v1_8_1:
	-y1_7_1 + x1_8_1 + y1_8_1 + z1_8_1 + y0_8_1 + z0_8_1 + x1_8_0 + y1_8_0 - 1.0 = 
	0;
subject to v2_8_1:
	-y2_7_1 - x1_8_1 + x2_8_1 + y2_8_1 + z2_8_1 + x2_8_0 + y2_8_0 - 1.0 = 0;
subject to v3_8_1:
	-y3_7_1 - x2_8_1 + x3_8_1 + y3_8_1 + z3_8_1 + x3_8_0 + y3_8_0 - 1.0 = 0;
subject to v4_8_1:
	-y4_7_1 - x3_8_1 + x4_8_1 + y4_8_1 + z4_8_1 + x4_8_0 + y4_8_0 - 1.0 = 0;
subject to v5_8_1:
	-y5_7_1 - x4_8_1 + x5_8_1 + y5_8_1 + z5_8_1 + x5_8_0 + y5_8_0 - 1.0 = 0;
subject to v6_8_1:
	-y6_7_1 - x5_8_1 + x6_8_1 + y6_8_1 + z6_8_1 + x6_8_0 + y6_8_0 - 1.0 = 0;
subject to v7_8_1:
	-y7_7_1 - x6_8_1 + x7_8_1 + y7_8_1 + z7_8_1 + x7_8_0 + y7_8_0 - 1.0 = 0;
subject to v8_8_1:
	-y8_7_1 - x7_8_1 + x8_8_1 + y8_8_1 + z8_8_1 + x8_8_0 + y8_8_0 - 1.0 = 0;
subject to v9_8_1:
	-y9_7_1 - x8_8_1 + x9_8_1 + y9_8_1 + z9_8_1 + x9_8_0 + y9_8_0 - 1.0 = 0;
subject to v10_8_1:
	-x9_8_1 - y10_7_1 + y10_8_1 + z10_8_1 + y11_8_1 + z11_8_1 + x10_8_0 + y10_8_0 - 
	1.0 = 0;
subject to v1_9_1:
	-y1_8_1 + x1_9_1 + y1_9_1 + z1_9_1 + y0_9_1 + z0_9_1 + x1_9_0 + y1_9_0 - 1.0 = 
	0;
subject to v2_9_1:
	-y2_8_1 - x1_9_1 + x2_9_1 + y2_9_1 + z2_9_1 + x2_9_0 + y2_9_0 - 1.0 = 0;
subject to v3_9_1:
	-y3_8_1 - x2_9_1 + x3_9_1 + y3_9_1 + z3_9_1 + x3_9_0 + y3_9_0 - 1.0 = 0;
subject to v4_9_1:
	-y4_8_1 - x3_9_1 + x4_9_1 + y4_9_1 + z4_9_1 + x4_9_0 + y4_9_0 - 1.0 = 0;
subject to v5_9_1:
	-y5_8_1 - x4_9_1 + x5_9_1 + y5_9_1 + z5_9_1 + x5_9_0 + y5_9_0 - 1.0 = 0;
subject to v6_9_1:
	-y6_8_1 - x5_9_1 + x6_9_1 + y6_9_1 + z6_9_1 + x6_9_0 + y6_9_0 - 1.0 = 0;
subject to v7_9_1:
	-y7_8_1 - x6_9_1 + x7_9_1 + y7_9_1 + z7_9_1 + x7_9_0 + y7_9_0 - 1.0 = 0;
subject to v8_9_1:
	-y8_8_1 - x7_9_1 + x8_9_1 + y8_9_1 + z8_9_1 + x8_9_0 + y8_9_0 - 1.0 = 0;
subject to v9_9_1:
	-y9_8_1 - x8_9_1 + x9_9_1 + y9_9_1 + z9_9_1 + x9_9_0 + y9_9_0 - 1.0 = 0;
subject to v10_9_1:
	-x9_9_1 - y10_8_1 + y10_9_1 + z10_9_1 + y11_9_1 + z11_9_1 + x10_9_0 + y10_9_0 - 
	1.0 = 0;
subject to v1_10_1:
	-y1_9_1 + x1_10_1 + z1_10_1 + y0_10_1 + z0_10_1 + x1_11_1 + z1_11_1 + x1_10_0 + 
	y1_10_0 - 1.0 = 0;
subject to v2_10_1:
	-y2_9_1 - x1_10_1 + x2_10_1 + z2_10_1 + x2_11_1 + z2_11_1 + x2_10_0 + y2_10_0 - 
	1.0 = 0;
subject to v3_10_1:
	-y3_9_1 - x2_10_1 + x3_10_1 + z3_10_1 + x3_11_1 + z3_11_1 + x3_10_0 + y3_10_0 - 
	1.0 = 0;
subject to v4_10_1:
	-y4_9_1 - x3_10_1 + x4_10_1 + z4_10_1 + x4_11_1 + z4_11_1 + x4_10_0 + y4_10_0 - 
	1.0 = 0;
subject to v5_10_1:
	-y5_9_1 - x4_10_1 + x5_10_1 + z5_10_1 + x5_11_1 + z5_11_1 + x5_10_0 + y5_10_0 - 
	1.0 = 0;
subject to v6_10_1:
	-y6_9_1 - x5_10_1 + x6_10_1 + z6_10_1 + x6_11_1 + z6_11_1 + x6_10_0 + y6_10_0 - 
	1.0 = 0;
subject to v7_10_1:
	-y7_9_1 - x6_10_1 + x7_10_1 + z7_10_1 + x7_11_1 + z7_11_1 + x7_10_0 + y7_10_0 - 
	1.0 = 0;
subject to v8_10_1:
	-y8_9_1 - x7_10_1 + x8_10_1 + z8_10_1 + x8_11_1 + z8_11_1 + x8_10_0 + y8_10_0 - 
	1.0 = 0;
subject to v9_10_1:
	-y9_9_1 - x8_10_1 + x9_10_1 + z9_10_1 + x9_11_1 + z9_11_1 + x9_10_0 + y9_10_0 - 
	1.0 = 0;
subject to v10_10_1:
	-y10_9_1 - x9_10_1 + y11_10_1 + z11_10_1 + x10_11_1 + z10_11_1 + x10_10_0 + 
	y10_10_0 - 1.0 = 0;
subject to v1_1_2:
	-z1_1_1 + x1_1_2 + y1_1_2 + z1_1_2 + y0_1_2 + z0_1_2 + x1_0_2 + z1_0_2 - 1.0 = 
	0;
subject to v2_1_2:
	-z2_1_1 - x1_1_2 + x2_1_2 + y2_1_2 + z2_1_2 + x2_0_2 + z2_0_2 - 1.0 = 0;
subject to v3_1_2:
	-z3_1_1 - x2_1_2 + x3_1_2 + y3_1_2 + z3_1_2 + x3_0_2 + z3_0_2 - 1.0 = 0;
subject to v4_1_2:
	-z4_1_1 - x3_1_2 + x4_1_2 + y4_1_2 + z4_1_2 + x4_0_2 + z4_0_2 - 1.0 = 0;
subject to v5_1_2:
	-z5_1_1 - x4_1_2 + x5_1_2 + y5_1_2 + z5_1_2 + x5_0_2 + z5_0_2 - 1.0 = 0;
subject to v6_1_2:
	-z6_1_1 - x5_1_2 + x6_1_2 + y6_1_2 + z6_1_2 + x6_0_2 + z6_0_2 - 1.0 = 0;
subject to v7_1_2:
	-z7_1_1 - x6_1_2 + x7_1_2 + y7_1_2 + z7_1_2 + x7_0_2 + z7_0_2 - 1.0 = 0;
subject to v8_1_2:
	-z8_1_1 - x7_1_2 + x8_1_2 + y8_1_2 + z8_1_2 + x8_0_2 + z8_0_2 - 1.0 = 0;
subject to v9_1_2:
	-z9_1_1 - x8_1_2 + x9_1_2 + y9_1_2 + z9_1_2 + x9_0_2 + z9_0_2 - 1.0 = 0;
subject to v10_1_2:
	-x9_1_2 - z10_1_1 + y10_1_2 + z10_1_2 + y11_1_2 + z11_1_2 + x10_0_2 + z10_0_2 - 
	1.0 = 0;
subject to v1_2_2:
	-z1_2_1 - y1_1_2 + x1_2_2 + y1_2_2 + z1_2_2 + y0_2_2 + z0_2_2 - 1.0 = 0;
subject to v2_2_2:
	-z2_2_1 - y2_1_2 - x1_2_2 + x2_2_2 + y2_2_2 + z2_2_2 - 1.0 = 0;
subject to v3_2_2:
	-z3_2_1 - y3_1_2 - x2_2_2 + x3_2_2 + y3_2_2 + z3_2_2 - 1.0 = 0;
subject to v4_2_2:
	-z4_2_1 - y4_1_2 - x3_2_2 + x4_2_2 + y4_2_2 + z4_2_2 - 1.0 = 0;
subject to v5_2_2:
	-z5_2_1 - y5_1_2 - x4_2_2 + x5_2_2 + y5_2_2 + z5_2_2 - 1.0 = 0;
subject to v6_2_2:
	-z6_2_1 - y6_1_2 - x5_2_2 + x6_2_2 + y6_2_2 + z6_2_2 - 1.0 = 0;
subject to v7_2_2:
	-z7_2_1 - y7_1_2 - x6_2_2 + x7_2_2 + y7_2_2 + z7_2_2 - 1.0 = 0;
subject to v8_2_2:
	-z8_2_1 - y8_1_2 - x7_2_2 + x8_2_2 + y8_2_2 + z8_2_2 - 1.0 = 0;
subject to v9_2_2:
	-z9_2_1 - y9_1_2 - x8_2_2 + x9_2_2 + y9_2_2 + z9_2_2 - 1.0 = 0;
subject to v10_2_2:
	-x9_2_2 - z10_2_1 - y10_1_2 + y10_2_2 + z10_2_2 + y11_2_2 + z11_2_2 - 1.0 = 0;
subject to v1_3_2:
	-z1_3_1 - y1_2_2 + x1_3_2 + y1_3_2 + z1_3_2 + y0_3_2 + z0_3_2 - 1.0 = 0;
subject to v2_3_2:
	-z2_3_1 - y2_2_2 - x1_3_2 + x2_3_2 + y2_3_2 + z2_3_2 - 1.0 = 0;
subject to v3_3_2:
	-z3_3_1 - y3_2_2 - x2_3_2 + x3_3_2 + y3_3_2 + z3_3_2 - 1.0 = 0;
subject to v4_3_2:
	-z4_3_1 - y4_2_2 - x3_3_2 + x4_3_2 + y4_3_2 + z4_3_2 - 1.0 = 0;
subject to v5_3_2:
	-z5_3_1 - y5_2_2 - x4_3_2 + x5_3_2 + y5_3_2 + z5_3_2 - 1.0 = 0;
subject to v6_3_2:
	-z6_3_1 - y6_2_2 - x5_3_2 + x6_3_2 + y6_3_2 + z6_3_2 - 1.0 = 0;
subject to v7_3_2:
	-z7_3_1 - y7_2_2 - x6_3_2 + x7_3_2 + y7_3_2 + z7_3_2 - 1.0 = 0;
subject to v8_3_2:
	-z8_3_1 - y8_2_2 - x7_3_2 + x8_3_2 + y8_3_2 + z8_3_2 - 1.0 = 0;
subject to v9_3_2:
	-z9_3_1 - y9_2_2 - x8_3_2 + x9_3_2 + y9_3_2 + z9_3_2 - 1.0 = 0;
subject to v10_3_2:
	-x9_3_2 - z10_3_1 - y10_2_2 + y10_3_2 + z10_3_2 + y11_3_2 + z11_3_2 - 1.0 = 0;
subject to v1_4_2:
	-z1_4_1 - y1_3_2 + x1_4_2 + y1_4_2 + z1_4_2 + y0_4_2 + z0_4_2 - 1.0 = 0;
subject to v2_4_2:
	-z2_4_1 - y2_3_2 - x1_4_2 + x2_4_2 + y2_4_2 + z2_4_2 - 1.0 = 0;
subject to v3_4_2:
	-z3_4_1 - y3_3_2 - x2_4_2 + x3_4_2 + y3_4_2 + z3_4_2 - 1.0 = 0;
subject to v4_4_2:
	-z4_4_1 - y4_3_2 - x3_4_2 + x4_4_2 + y4_4_2 + z4_4_2 - 1.0 = 0;
subject to v5_4_2:
	-z5_4_1 - y5_3_2 - x4_4_2 + x5_4_2 + y5_4_2 + z5_4_2 - 1.0 = 0;
subject to v6_4_2:
	-z6_4_1 - y6_3_2 - x5_4_2 + x6_4_2 + y6_4_2 + z6_4_2 - 1.0 = 0;
subject to v7_4_2:
	-z7_4_1 - y7_3_2 - x6_4_2 + x7_4_2 + y7_4_2 + z7_4_2 - 1.0 = 0;
subject to v8_4_2:
	-z8_4_1 - y8_3_2 - x7_4_2 + x8_4_2 + y8_4_2 + z8_4_2 - 1.0 = 0;
subject to v9_4_2:
	-z9_4_1 - y9_3_2 - x8_4_2 + x9_4_2 + y9_4_2 + z9_4_2 - 1.0 = 0;
subject to v10_4_2:
	-x9_4_2 - z10_4_1 - y10_3_2 + y10_4_2 + z10_4_2 + y11_4_2 + z11_4_2 - 1.0 = 0;
subject to v1_5_2:
	-z1_5_1 - y1_4_2 + x1_5_2 + y1_5_2 + z1_5_2 + y0_5_2 + z0_5_2 - 1.0 = 0;
subject to v2_5_2:
	-z2_5_1 - y2_4_2 - x1_5_2 + x2_5_2 + y2_5_2 + z2_5_2 - 1.0 = 0;
subject to v3_5_2:
	-z3_5_1 - y3_4_2 - x2_5_2 + x3_5_2 + y3_5_2 + z3_5_2 - 1.0 = 0;
subject to v4_5_2:
	-z4_5_1 - y4_4_2 - x3_5_2 + x4_5_2 + y4_5_2 + z4_5_2 - 1.0 = 0;
subject to v5_5_2:
	-z5_5_1 - y5_4_2 - x4_5_2 + x5_5_2 + y5_5_2 + z5_5_2 - 1.0 = 0;
subject to v6_5_2:
	-z6_5_1 - y6_4_2 - x5_5_2 + x6_5_2 + y6_5_2 + z6_5_2 - 1.0 = 0;
subject to v7_5_2:
	-z7_5_1 - y7_4_2 - x6_5_2 + x7_5_2 + y7_5_2 + z7_5_2 - 1.0 = 0;
subject to v8_5_2:
	-z8_5_1 - y8_4_2 - x7_5_2 + x8_5_2 + y8_5_2 + z8_5_2 - 1.0 = 0;
subject to v9_5_2:
	-z9_5_1 - y9_4_2 - x8_5_2 + x9_5_2 + y9_5_2 + z9_5_2 - 1.0 = 0;
subject to v10_5_2:
	-x9_5_2 - z10_5_1 - y10_4_2 + y10_5_2 + z10_5_2 + y11_5_2 + z11_5_2 - 1.0 = 0;
subject to v1_6_2:
	-z1_6_1 - y1_5_2 + x1_6_2 + y1_6_2 + z1_6_2 + y0_6_2 + z0_6_2 - 1.0 = 0;
subject to v2_6_2:
	-z2_6_1 - y2_5_2 - x1_6_2 + x2_6_2 + y2_6_2 + z2_6_2 - 1.0 = 0;
subject to v3_6_2:
	-z3_6_1 - y3_5_2 - x2_6_2 + x3_6_2 + y3_6_2 + z3_6_2 - 1.0 = 0;
subject to v4_6_2:
	-z4_6_1 - y4_5_2 - x3_6_2 + x4_6_2 + y4_6_2 + z4_6_2 - 1.0 = 0;
subject to v5_6_2:
	-z5_6_1 - y5_5_2 - x4_6_2 + x5_6_2 + y5_6_2 + z5_6_2 - 1.0 = 0;
subject to v6_6_2:
	-z6_6_1 - y6_5_2 - x5_6_2 + x6_6_2 + y6_6_2 + z6_6_2 - 1.0 = 0;
subject to v7_6_2:
	-z7_6_1 - y7_5_2 - x6_6_2 + x7_6_2 + y7_6_2 + z7_6_2 - 1.0 = 0;
subject to v8_6_2:
	-z8_6_1 - y8_5_2 - x7_6_2 + x8_6_2 + y8_6_2 + z8_6_2 - 1.0 = 0;
subject to v9_6_2:
	-z9_6_1 - y9_5_2 - x8_6_2 + x9_6_2 + y9_6_2 + z9_6_2 - 1.0 = 0;
subject to v10_6_2:
	-x9_6_2 - z10_6_1 - y10_5_2 + y10_6_2 + z10_6_2 + y11_6_2 + z11_6_2 - 1.0 = 0;
subject to v1_7_2:
	-z1_7_1 - y1_6_2 + x1_7_2 + y1_7_2 + z1_7_2 + y0_7_2 + z0_7_2 - 1.0 = 0;
subject to v2_7_2:
	-z2_7_1 - y2_6_2 - x1_7_2 + x2_7_2 + y2_7_2 + z2_7_2 - 1.0 = 0;
subject to v3_7_2:
	-z3_7_1 - y3_6_2 - x2_7_2 + x3_7_2 + y3_7_2 + z3_7_2 - 1.0 = 0;
subject to v4_7_2:
	-z4_7_1 - y4_6_2 - x3_7_2 + x4_7_2 + y4_7_2 + z4_7_2 - 1.0 = 0;
subject to v5_7_2:
	-z5_7_1 - y5_6_2 - x4_7_2 + x5_7_2 + y5_7_2 + z5_7_2 - 1.0 = 0;
subject to v6_7_2:
	-z6_7_1 - y6_6_2 - x5_7_2 + x6_7_2 + y6_7_2 + z6_7_2 - 1.0 = 0;
subject to v7_7_2:
	-z7_7_1 - y7_6_2 - x6_7_2 + x7_7_2 + y7_7_2 + z7_7_2 - 1.0 = 0;
subject to v8_7_2:
	-z8_7_1 - y8_6_2 - x7_7_2 + x8_7_2 + y8_7_2 + z8_7_2 - 1.0 = 0;
subject to v9_7_2:
	-z9_7_1 - y9_6_2 - x8_7_2 + x9_7_2 + y9_7_2 + z9_7_2 - 1.0 = 0;
subject to v10_7_2:
	-x9_7_2 - z10_7_1 - y10_6_2 + y10_7_2 + z10_7_2 + y11_7_2 + z11_7_2 - 1.0 = 0;
subject to v1_8_2:
	-z1_8_1 - y1_7_2 + x1_8_2 + y1_8_2 + z1_8_2 + y0_8_2 + z0_8_2 - 1.0 = 0;
subject to v2_8_2:
	-z2_8_1 - y2_7_2 - x1_8_2 + x2_8_2 + y2_8_2 + z2_8_2 - 1.0 = 0;
subject to v3_8_2:
	-z3_8_1 - y3_7_2 - x2_8_2 + x3_8_2 + y3_8_2 + z3_8_2 - 1.0 = 0;
subject to v4_8_2:
	-z4_8_1 - y4_7_2 - x3_8_2 + x4_8_2 + y4_8_2 + z4_8_2 - 1.0 = 0;
subject to v5_8_2:
	-z5_8_1 - y5_7_2 - x4_8_2 + x5_8_2 + y5_8_2 + z5_8_2 - 1.0 = 0;
subject to v6_8_2:
	-z6_8_1 - y6_7_2 - x5_8_2 + x6_8_2 + y6_8_2 + z6_8_2 - 1.0 = 0;
subject to v7_8_2:
	-z7_8_1 - y7_7_2 - x6_8_2 + x7_8_2 + y7_8_2 + z7_8_2 - 1.0 = 0;
subject to v8_8_2:
	-z8_8_1 - y8_7_2 - x7_8_2 + x8_8_2 + y8_8_2 + z8_8_2 - 1.0 = 0;
subject to v9_8_2:
	-z9_8_1 - y9_7_2 - x8_8_2 + x9_8_2 + y9_8_2 + z9_8_2 - 1.0 = 0;
subject to v10_8_2:
	-x9_8_2 - z10_8_1 - y10_7_2 + y10_8_2 + z10_8_2 + y11_8_2 + z11_8_2 - 1.0 = 0;
subject to v1_9_2:
	-z1_9_1 - y1_8_2 + x1_9_2 + y1_9_2 + z1_9_2 + y0_9_2 + z0_9_2 - 1.0 = 0;
subject to v2_9_2:
	-z2_9_1 - y2_8_2 - x1_9_2 + x2_9_2 + y2_9_2 + z2_9_2 - 1.0 = 0;
subject to v3_9_2:
	-z3_9_1 - y3_8_2 - x2_9_2 + x3_9_2 + y3_9_2 + z3_9_2 - 1.0 = 0;
subject to v4_9_2:
	-z4_9_1 - y4_8_2 - x3_9_2 + x4_9_2 + y4_9_2 + z4_9_2 - 1.0 = 0;
subject to v5_9_2:
	-z5_9_1 - y5_8_2 - x4_9_2 + x5_9_2 + y5_9_2 + z5_9_2 - 1.0 = 0;
subject to v6_9_2:
	-z6_9_1 - y6_8_2 - x5_9_2 + x6_9_2 + y6_9_2 + z6_9_2 - 1.0 = 0;
subject to v7_9_2:
	-z7_9_1 - y7_8_2 - x6_9_2 + x7_9_2 + y7_9_2 + z7_9_2 - 1.0 = 0;
subject to v8_9_2:
	-z8_9_1 - y8_8_2 - x7_9_2 + x8_9_2 + y8_9_2 + z8_9_2 - 1.0 = 0;
subject to v9_9_2:
	-z9_9_1 - y9_8_2 - x8_9_2 + x9_9_2 + y9_9_2 + z9_9_2 - 1.0 = 0;
subject to v10_9_2:
	-x9_9_2 - z10_9_1 - y10_8_2 + y10_9_2 + z10_9_2 + y11_9_2 + z11_9_2 - 1.0 = 0;
subject to v1_10_2:
	-y1_9_2 - z1_10_1 + x1_10_2 + z1_10_2 + y0_10_2 + z0_10_2 + x1_11_2 + z1_11_2 - 
	1.0 = 0;
subject to v2_10_2:
	-y2_9_2 - z2_10_1 - x1_10_2 + x2_10_2 + z2_10_2 + x2_11_2 + z2_11_2 - 1.0 = 0;
subject to v3_10_2:
	-y3_9_2 - z3_10_1 - x2_10_2 + x3_10_2 + z3_10_2 + x3_11_2 + z3_11_2 - 1.0 = 0;
subject to v4_10_2:
	-y4_9_2 - z4_10_1 - x3_10_2 + x4_10_2 + z4_10_2 + x4_11_2 + z4_11_2 - 1.0 = 0;
subject to v5_10_2:
	-y5_9_2 - z5_10_1 - x4_10_2 + x5_10_2 + z5_10_2 + x5_11_2 + z5_11_2 - 1.0 = 0;
subject to v6_10_2:
	-y6_9_2 - z6_10_1 - x5_10_2 + x6_10_2 + z6_10_2 + x6_11_2 + z6_11_2 - 1.0 = 0;
subject to v7_10_2:
	-y7_9_2 - z7_10_1 - x6_10_2 + x7_10_2 + z7_10_2 + x7_11_2 + z7_11_2 - 1.0 = 0;
subject to v8_10_2:
	-y8_9_2 - z8_10_1 - x7_10_2 + x8_10_2 + z8_10_2 + x8_11_2 + z8_11_2 - 1.0 = 0;
subject to v9_10_2:
	-y9_9_2 - z9_10_1 - x8_10_2 + x9_10_2 + z9_10_2 + x9_11_2 + z9_11_2 - 1.0 = 0;
subject to v10_10_2:
	-y10_9_2 - x9_10_2 + y11_10_2 + z11_10_2 + x10_11_2 + z10_11_2 - 1.0 = 0;
subject to v1_1_3:
	-z1_1_2 + x1_1_3 + y1_1_3 + z1_1_3 + y0_1_3 + z0_1_3 + x1_0_3 + z1_0_3 - 1.0 = 
	0;
subject to v2_1_3:
	-z2_1_2 - x1_1_3 + x2_1_3 + y2_1_3 + z2_1_3 + x2_0_3 + z2_0_3 - 1.0 = 0;
subject to v3_1_3:
	-z3_1_2 - x2_1_3 + x3_1_3 + y3_1_3 + z3_1_3 + x3_0_3 + z3_0_3 - 1.0 = 0;
subject to v4_1_3:
	-z4_1_2 - x3_1_3 + x4_1_3 + y4_1_3 + z4_1_3 + x4_0_3 + z4_0_3 - 1.0 = 0;
subject to v5_1_3:
	-z5_1_2 - x4_1_3 + x5_1_3 + y5_1_3 + z5_1_3 + x5_0_3 + z5_0_3 - 1.0 = 0;
subject to v6_1_3:
	-z6_1_2 - x5_1_3 + x6_1_3 + y6_1_3 + z6_1_3 + x6_0_3 + z6_0_3 - 1.0 = 0;
subject to v7_1_3:
	-z7_1_2 - x6_1_3 + x7_1_3 + y7_1_3 + z7_1_3 + x7_0_3 + z7_0_3 - 1.0 = 0;
subject to v8_1_3:
	-z8_1_2 - x7_1_3 + x8_1_3 + y8_1_3 + z8_1_3 + x8_0_3 + z8_0_3 - 1.0 = 0;
subject to v9_1_3:
	-z9_1_2 - x8_1_3 + x9_1_3 + y9_1_3 + z9_1_3 + x9_0_3 + z9_0_3 - 1.0 = 0;
subject to v10_1_3:
	-x9_1_3 - z10_1_2 + y10_1_3 + z10_1_3 + y11_1_3 + z11_1_3 + x10_0_3 + z10_0_3 - 
	1.0 = 0;
subject to v1_2_3:
	-z1_2_2 - y1_1_3 + x1_2_3 + y1_2_3 + z1_2_3 + y0_2_3 + z0_2_3 - 1.0 = 0;
subject to v2_2_3:
	-z2_2_2 - y2_1_3 - x1_2_3 + x2_2_3 + y2_2_3 + z2_2_3 - 1.0 = 0;
subject to v3_2_3:
	-z3_2_2 - y3_1_3 - x2_2_3 + x3_2_3 + y3_2_3 + z3_2_3 - 1.0 = 0;
subject to v4_2_3:
	-z4_2_2 - y4_1_3 - x3_2_3 + x4_2_3 + y4_2_3 + z4_2_3 - 1.0 = 0;
subject to v5_2_3:
	-z5_2_2 - y5_1_3 - x4_2_3 + x5_2_3 + y5_2_3 + z5_2_3 - 1.0 = 0;
subject to v6_2_3:
	-z6_2_2 - y6_1_3 - x5_2_3 + x6_2_3 + y6_2_3 + z6_2_3 - 1.0 = 0;
subject to v7_2_3:
	-z7_2_2 - y7_1_3 - x6_2_3 + x7_2_3 + y7_2_3 + z7_2_3 - 1.0 = 0;
subject to v8_2_3:
	-z8_2_2 - y8_1_3 - x7_2_3 + x8_2_3 + y8_2_3 + z8_2_3 - 1.0 = 0;
subject to v9_2_3:
	-z9_2_2 - y9_1_3 - x8_2_3 + x9_2_3 + y9_2_3 + z9_2_3 - 1.0 = 0;
subject to v10_2_3:
	-x9_2_3 - z10_2_2 - y10_1_3 + y10_2_3 + z10_2_3 + y11_2_3 + z11_2_3 - 1.0 = 0;
subject to v1_3_3:
	-z1_3_2 - y1_2_3 + x1_3_3 + y1_3_3 + z1_3_3 + y0_3_3 + z0_3_3 - 1.0 = 0;
subject to v2_3_3:
	-z2_3_2 - y2_2_3 - x1_3_3 + x2_3_3 + y2_3_3 + z2_3_3 - 1.0 = 0;
subject to v3_3_3:
	-z3_3_2 - y3_2_3 - x2_3_3 + x3_3_3 + y3_3_3 + z3_3_3 - 1.0 = 0;
subject to v4_3_3:
	-z4_3_2 - y4_2_3 - x3_3_3 + x4_3_3 + y4_3_3 + z4_3_3 - 1.0 = 0;
subject to v5_3_3:
	-z5_3_2 - y5_2_3 - x4_3_3 + x5_3_3 + y5_3_3 + z5_3_3 - 1.0 = 0;
subject to v6_3_3:
	-z6_3_2 - y6_2_3 - x5_3_3 + x6_3_3 + y6_3_3 + z6_3_3 - 1.0 = 0;
subject to v7_3_3:
	-z7_3_2 - y7_2_3 - x6_3_3 + x7_3_3 + y7_3_3 + z7_3_3 - 1.0 = 0;
subject to v8_3_3:
	-z8_3_2 - y8_2_3 - x7_3_3 + x8_3_3 + y8_3_3 + z8_3_3 - 1.0 = 0;
subject to v9_3_3:
	-z9_3_2 - y9_2_3 - x8_3_3 + x9_3_3 + y9_3_3 + z9_3_3 - 1.0 = 0;
subject to v10_3_3:
	-x9_3_3 - z10_3_2 - y10_2_3 + y10_3_3 + z10_3_3 + y11_3_3 + z11_3_3 - 1.0 = 0;
subject to v1_4_3:
	-z1_4_2 - y1_3_3 + x1_4_3 + y1_4_3 + z1_4_3 + y0_4_3 + z0_4_3 - 1.0 = 0;
subject to v2_4_3:
	-z2_4_2 - y2_3_3 - x1_4_3 + x2_4_3 + y2_4_3 + z2_4_3 - 1.0 = 0;
subject to v3_4_3:
	-z3_4_2 - y3_3_3 - x2_4_3 + x3_4_3 + y3_4_3 + z3_4_3 - 1.0 = 0;
subject to v4_4_3:
	-z4_4_2 - y4_3_3 - x3_4_3 + x4_4_3 + y4_4_3 + z4_4_3 - 1.0 = 0;
subject to v5_4_3:
	-z5_4_2 - y5_3_3 - x4_4_3 + x5_4_3 + y5_4_3 + z5_4_3 - 1.0 = 0;
subject to v6_4_3:
	-z6_4_2 - y6_3_3 - x5_4_3 + x6_4_3 + y6_4_3 + z6_4_3 - 1.0 = 0;
subject to v7_4_3:
	-z7_4_2 - y7_3_3 - x6_4_3 + x7_4_3 + y7_4_3 + z7_4_3 - 1.0 = 0;
subject to v8_4_3:
	-z8_4_2 - y8_3_3 - x7_4_3 + x8_4_3 + y8_4_3 + z8_4_3 - 1.0 = 0;
subject to v9_4_3:
	-z9_4_2 - y9_3_3 - x8_4_3 + x9_4_3 + y9_4_3 + z9_4_3 - 1.0 = 0;
subject to v10_4_3:
	-x9_4_3 - z10_4_2 - y10_3_3 + y10_4_3 + z10_4_3 + y11_4_3 + z11_4_3 - 1.0 = 0;
subject to v1_5_3:
	-z1_5_2 - y1_4_3 + x1_5_3 + y1_5_3 + z1_5_3 + y0_5_3 + z0_5_3 - 1.0 = 0;
subject to v2_5_3:
	-z2_5_2 - y2_4_3 - x1_5_3 + x2_5_3 + y2_5_3 + z2_5_3 - 1.0 = 0;
subject to v3_5_3:
	-z3_5_2 - y3_4_3 - x2_5_3 + x3_5_3 + y3_5_3 + z3_5_3 - 1.0 = 0;
subject to v4_5_3:
	-z4_5_2 - y4_4_3 - x3_5_3 + x4_5_3 + y4_5_3 + z4_5_3 - 1.0 = 0;
subject to v5_5_3:
	-z5_5_2 - y5_4_3 - x4_5_3 + x5_5_3 + y5_5_3 + z5_5_3 - 1.0 = 0;
subject to v6_5_3:
	-z6_5_2 - y6_4_3 - x5_5_3 + x6_5_3 + y6_5_3 + z6_5_3 - 1.0 = 0;
subject to v7_5_3:
	-z7_5_2 - y7_4_3 - x6_5_3 + x7_5_3 + y7_5_3 + z7_5_3 - 1.0 = 0;
subject to v8_5_3:
	-z8_5_2 - y8_4_3 - x7_5_3 + x8_5_3 + y8_5_3 + z8_5_3 - 1.0 = 0;
subject to v9_5_3:
	-z9_5_2 - y9_4_3 - x8_5_3 + x9_5_3 + y9_5_3 + z9_5_3 - 1.0 = 0;
subject to v10_5_3:
	-x9_5_3 - z10_5_2 - y10_4_3 + y10_5_3 + z10_5_3 + y11_5_3 + z11_5_3 - 1.0 = 0;
subject to v1_6_3:
	-z1_6_2 - y1_5_3 + x1_6_3 + y1_6_3 + z1_6_3 + y0_6_3 + z0_6_3 - 1.0 = 0;
subject to v2_6_3:
	-z2_6_2 - y2_5_3 - x1_6_3 + x2_6_3 + y2_6_3 + z2_6_3 - 1.0 = 0;
subject to v3_6_3:
	-z3_6_2 - y3_5_3 - x2_6_3 + x3_6_3 + y3_6_3 + z3_6_3 - 1.0 = 0;
subject to v4_6_3:
	-z4_6_2 - y4_5_3 - x3_6_3 + x4_6_3 + y4_6_3 + z4_6_3 - 1.0 = 0;
subject to v5_6_3:
	-z5_6_2 - y5_5_3 - x4_6_3 + x5_6_3 + y5_6_3 + z5_6_3 - 1.0 = 0;
subject to v6_6_3:
	-z6_6_2 - y6_5_3 - x5_6_3 + x6_6_3 + y6_6_3 + z6_6_3 - 1.0 = 0;
subject to v7_6_3:
	-z7_6_2 - y7_5_3 - x6_6_3 + x7_6_3 + y7_6_3 + z7_6_3 - 1.0 = 0;
subject to v8_6_3:
	-z8_6_2 - y8_5_3 - x7_6_3 + x8_6_3 + y8_6_3 + z8_6_3 - 1.0 = 0;
subject to v9_6_3:
	-z9_6_2 - y9_5_3 - x8_6_3 + x9_6_3 + y9_6_3 + z9_6_3 - 1.0 = 0;
subject to v10_6_3:
	-x9_6_3 - z10_6_2 - y10_5_3 + y10_6_3 + z10_6_3 + y11_6_3 + z11_6_3 - 1.0 = 0;
subject to v1_7_3:
	-z1_7_2 - y1_6_3 + x1_7_3 + y1_7_3 + z1_7_3 + y0_7_3 + z0_7_3 - 1.0 = 0;
subject to v2_7_3:
	-z2_7_2 - y2_6_3 - x1_7_3 + x2_7_3 + y2_7_3 + z2_7_3 - 1.0 = 0;
subject to v3_7_3:
	-z3_7_2 - y3_6_3 - x2_7_3 + x3_7_3 + y3_7_3 + z3_7_3 - 1.0 = 0;
subject to v4_7_3:
	-z4_7_2 - y4_6_3 - x3_7_3 + x4_7_3 + y4_7_3 + z4_7_3 - 1.0 = 0;
subject to v5_7_3:
	-z5_7_2 - y5_6_3 - x4_7_3 + x5_7_3 + y5_7_3 + z5_7_3 - 1.0 = 0;
subject to v6_7_3:
	-z6_7_2 - y6_6_3 - x5_7_3 + x6_7_3 + y6_7_3 + z6_7_3 - 1.0 = 0;
subject to v7_7_3:
	-z7_7_2 - y7_6_3 - x6_7_3 + x7_7_3 + y7_7_3 + z7_7_3 - 1.0 = 0;
subject to v8_7_3:
	-z8_7_2 - y8_6_3 - x7_7_3 + x8_7_3 + y8_7_3 + z8_7_3 - 1.0 = 0;
subject to v9_7_3:
	-z9_7_2 - y9_6_3 - x8_7_3 + x9_7_3 + y9_7_3 + z9_7_3 - 1.0 = 0;
subject to v10_7_3:
	-x9_7_3 - z10_7_2 - y10_6_3 + y10_7_3 + z10_7_3 + y11_7_3 + z11_7_3 - 1.0 = 0;
subject to v1_8_3:
	-z1_8_2 - y1_7_3 + x1_8_3 + y1_8_3 + z1_8_3 + y0_8_3 + z0_8_3 - 1.0 = 0;
subject to v2_8_3:
	-z2_8_2 - y2_7_3 - x1_8_3 + x2_8_3 + y2_8_3 + z2_8_3 - 1.0 = 0;
subject to v3_8_3:
	-z3_8_2 - y3_7_3 - x2_8_3 + x3_8_3 + y3_8_3 + z3_8_3 - 1.0 = 0;
subject to v4_8_3:
	-z4_8_2 - y4_7_3 - x3_8_3 + x4_8_3 + y4_8_3 + z4_8_3 - 1.0 = 0;
subject to v5_8_3:
	-z5_8_2 - y5_7_3 - x4_8_3 + x5_8_3 + y5_8_3 + z5_8_3 - 1.0 = 0;
subject to v6_8_3:
	-z6_8_2 - y6_7_3 - x5_8_3 + x6_8_3 + y6_8_3 + z6_8_3 - 1.0 = 0;
subject to v7_8_3:
	-z7_8_2 - y7_7_3 - x6_8_3 + x7_8_3 + y7_8_3 + z7_8_3 - 1.0 = 0;
subject to v8_8_3:
	-z8_8_2 - y8_7_3 - x7_8_3 + x8_8_3 + y8_8_3 + z8_8_3 - 1.0 = 0;
subject to v9_8_3:
	-z9_8_2 - y9_7_3 - x8_8_3 + x9_8_3 + y9_8_3 + z9_8_3 - 1.0 = 0;
subject to v10_8_3:
	-x9_8_3 - z10_8_2 - y10_7_3 + y10_8_3 + z10_8_3 + y11_8_3 + z11_8_3 - 1.0 = 0;
subject to v1_9_3:
	-z1_9_2 - y1_8_3 + x1_9_3 + y1_9_3 + z1_9_3 + y0_9_3 + z0_9_3 - 1.0 = 0;
subject to v2_9_3:
	-z2_9_2 - y2_8_3 - x1_9_3 + x2_9_3 + y2_9_3 + z2_9_3 - 1.0 = 0;
subject to v3_9_3:
	-z3_9_2 - y3_8_3 - x2_9_3 + x3_9_3 + y3_9_3 + z3_9_3 - 1.0 = 0;
subject to v4_9_3:
	-z4_9_2 - y4_8_3 - x3_9_3 + x4_9_3 + y4_9_3 + z4_9_3 - 1.0 = 0;
subject to v5_9_3:
	-z5_9_2 - y5_8_3 - x4_9_3 + x5_9_3 + y5_9_3 + z5_9_3 - 1.0 = 0;
subject to v6_9_3:
	-z6_9_2 - y6_8_3 - x5_9_3 + x6_9_3 + y6_9_3 + z6_9_3 - 1.0 = 0;
subject to v7_9_3:
	-z7_9_2 - y7_8_3 - x6_9_3 + x7_9_3 + y7_9_3 + z7_9_3 - 1.0 = 0;
subject to v8_9_3:
	-z8_9_2 - y8_8_3 - x7_9_3 + x8_9_3 + y8_9_3 + z8_9_3 - 1.0 = 0;
subject to v9_9_3:
	-z9_9_2 - y9_8_3 - x8_9_3 + x9_9_3 + y9_9_3 + z9_9_3 - 1.0 = 0;
subject to v10_9_3:
	-x9_9_3 - z10_9_2 - y10_8_3 + y10_9_3 + z10_9_3 + y11_9_3 + z11_9_3 - 1.0 = 0;
subject to v1_10_3:
	-y1_9_3 - z1_10_2 + x1_10_3 + z1_10_3 + y0_10_3 + z0_10_3 + x1_11_3 + z1_11_3 - 
	1.0 = 0;
subject to v2_10_3:
	-y2_9_3 - z2_10_2 - x1_10_3 + x2_10_3 + z2_10_3 + x2_11_3 + z2_11_3 - 1.0 = 0;
subject to v3_10_3:
	-y3_9_3 - z3_10_2 - x2_10_3 + x3_10_3 + z3_10_3 + x3_11_3 + z3_11_3 - 1.0 = 0;
subject to v4_10_3:
	-y4_9_3 - z4_10_2 - x3_10_3 + x4_10_3 + z4_10_3 + x4_11_3 + z4_11_3 - 1.0 = 0;
subject to v5_10_3:
	-y5_9_3 - z5_10_2 - x4_10_3 + x5_10_3 + z5_10_3 + x5_11_3 + z5_11_3 - 1.0 = 0;
subject to v6_10_3:
	-y6_9_3 - z6_10_2 - x5_10_3 + x6_10_3 + z6_10_3 + x6_11_3 + z6_11_3 - 1.0 = 0;
subject to v7_10_3:
	-y7_9_3 - z7_10_2 - x6_10_3 + x7_10_3 + z7_10_3 + x7_11_3 + z7_11_3 - 1.0 = 0;
subject to v8_10_3:
	-y8_9_3 - z8_10_2 - x7_10_3 + x8_10_3 + z8_10_3 + x8_11_3 + z8_11_3 - 1.0 = 0;
subject to v9_10_3:
	-y9_9_3 - z9_10_2 - x8_10_3 + x9_10_3 + z9_10_3 + x9_11_3 + z9_11_3 - 1.0 = 0;
subject to v10_10_3:
	-y10_9_3 - x9_10_3 + y11_10_3 + z11_10_3 + x10_11_3 + z10_11_3 - 1.0 = 0;
subject to v1_1_4:
	-z1_1_3 + x1_1_4 + y1_1_4 + z1_1_4 + y0_1_4 + z0_1_4 + x1_0_4 + z1_0_4 - 1.0 = 
	0;
subject to v2_1_4:
	-z2_1_3 - x1_1_4 + x2_1_4 + y2_1_4 + z2_1_4 + x2_0_4 + z2_0_4 - 1.0 = 0;
subject to v3_1_4:
	-z3_1_3 - x2_1_4 + x3_1_4 + y3_1_4 + z3_1_4 + x3_0_4 + z3_0_4 - 1.0 = 0;
subject to v4_1_4:
	-z4_1_3 - x3_1_4 + x4_1_4 + y4_1_4 + z4_1_4 + x4_0_4 + z4_0_4 - 1.0 = 0;
subject to v5_1_4:
	-z5_1_3 - x4_1_4 + x5_1_4 + y5_1_4 + z5_1_4 + x5_0_4 + z5_0_4 - 1.0 = 0;
subject to v6_1_4:
	-z6_1_3 - x5_1_4 + x6_1_4 + y6_1_4 + z6_1_4 + x6_0_4 + z6_0_4 - 1.0 = 0;
subject to v7_1_4:
	-z7_1_3 - x6_1_4 + x7_1_4 + y7_1_4 + z7_1_4 + x7_0_4 + z7_0_4 - 1.0 = 0;
subject to v8_1_4:
	-z8_1_3 - x7_1_4 + x8_1_4 + y8_1_4 + z8_1_4 + x8_0_4 + z8_0_4 - 1.0 = 0;
subject to v9_1_4:
	-z9_1_3 - x8_1_4 + x9_1_4 + y9_1_4 + z9_1_4 + x9_0_4 + z9_0_4 - 1.0 = 0;
subject to v10_1_4:
	-x9_1_4 - z10_1_3 + y10_1_4 + z10_1_4 + y11_1_4 + z11_1_4 + x10_0_4 + z10_0_4 - 
	1.0 = 0;
subject to v1_2_4:
	-z1_2_3 - y1_1_4 + x1_2_4 + y1_2_4 + z1_2_4 + y0_2_4 + z0_2_4 - 1.0 = 0;
subject to v2_2_4:
	-z2_2_3 - y2_1_4 - x1_2_4 + x2_2_4 + y2_2_4 + z2_2_4 - 1.0 = 0;
subject to v3_2_4:
	-z3_2_3 - y3_1_4 - x2_2_4 + x3_2_4 + y3_2_4 + z3_2_4 - 1.0 = 0;
subject to v4_2_4:
	-z4_2_3 - y4_1_4 - x3_2_4 + x4_2_4 + y4_2_4 + z4_2_4 - 1.0 = 0;
subject to v5_2_4:
	-z5_2_3 - y5_1_4 - x4_2_4 + x5_2_4 + y5_2_4 + z5_2_4 - 1.0 = 0;
subject to v6_2_4:
	-z6_2_3 - y6_1_4 - x5_2_4 + x6_2_4 + y6_2_4 + z6_2_4 - 1.0 = 0;
subject to v7_2_4:
	-z7_2_3 - y7_1_4 - x6_2_4 + x7_2_4 + y7_2_4 + z7_2_4 - 1.0 = 0;
subject to v8_2_4:
	-z8_2_3 - y8_1_4 - x7_2_4 + x8_2_4 + y8_2_4 + z8_2_4 - 1.0 = 0;
subject to v9_2_4:
	-z9_2_3 - y9_1_4 - x8_2_4 + x9_2_4 + y9_2_4 + z9_2_4 - 1.0 = 0;
subject to v10_2_4:
	-x9_2_4 - z10_2_3 - y10_1_4 + y10_2_4 + z10_2_4 + y11_2_4 + z11_2_4 - 1.0 = 0;
subject to v1_3_4:
	-z1_3_3 - y1_2_4 + x1_3_4 + y1_3_4 + z1_3_4 + y0_3_4 + z0_3_4 - 1.0 = 0;
subject to v2_3_4:
	-z2_3_3 - y2_2_4 - x1_3_4 + x2_3_4 + y2_3_4 + z2_3_4 - 1.0 = 0;
subject to v3_3_4:
	-z3_3_3 - y3_2_4 - x2_3_4 + x3_3_4 + y3_3_4 + z3_3_4 - 1.0 = 0;
subject to v4_3_4:
	-z4_3_3 - y4_2_4 - x3_3_4 + x4_3_4 + y4_3_4 + z4_3_4 - 1.0 = 0;
subject to v5_3_4:
	-z5_3_3 - y5_2_4 - x4_3_4 + x5_3_4 + y5_3_4 + z5_3_4 - 1.0 = 0;
subject to v6_3_4:
	-z6_3_3 - y6_2_4 - x5_3_4 + x6_3_4 + y6_3_4 + z6_3_4 - 1.0 = 0;
subject to v7_3_4:
	-z7_3_3 - y7_2_4 - x6_3_4 + x7_3_4 + y7_3_4 + z7_3_4 - 1.0 = 0;
subject to v8_3_4:
	-z8_3_3 - y8_2_4 - x7_3_4 + x8_3_4 + y8_3_4 + z8_3_4 - 1.0 = 0;
subject to v9_3_4:
	-z9_3_3 - y9_2_4 - x8_3_4 + x9_3_4 + y9_3_4 + z9_3_4 - 1.0 = 0;
subject to v10_3_4:
	-x9_3_4 - z10_3_3 - y10_2_4 + y10_3_4 + z10_3_4 + y11_3_4 + z11_3_4 - 1.0 = 0;
subject to v1_4_4:
	-z1_4_3 - y1_3_4 + x1_4_4 + y1_4_4 + z1_4_4 + y0_4_4 + z0_4_4 - 1.0 = 0;
subject to v2_4_4:
	-z2_4_3 - y2_3_4 - x1_4_4 + x2_4_4 + y2_4_4 + z2_4_4 - 1.0 = 0;
subject to v3_4_4:
	-z3_4_3 - y3_3_4 - x2_4_4 + x3_4_4 + y3_4_4 + z3_4_4 - 1.0 = 0;
subject to v4_4_4:
	-z4_4_3 - y4_3_4 - x3_4_4 + x4_4_4 + y4_4_4 + z4_4_4 - 1.0 = 0;
subject to v5_4_4:
	-z5_4_3 - y5_3_4 - x4_4_4 + x5_4_4 + y5_4_4 + z5_4_4 - 1.0 = 0;
subject to v6_4_4:
	-z6_4_3 - y6_3_4 - x5_4_4 + x6_4_4 + y6_4_4 + z6_4_4 - 1.0 = 0;
subject to v7_4_4:
	-z7_4_3 - y7_3_4 - x6_4_4 + x7_4_4 + y7_4_4 + z7_4_4 - 1.0 = 0;
subject to v8_4_4:
	-z8_4_3 - y8_3_4 - x7_4_4 + x8_4_4 + y8_4_4 + z8_4_4 - 1.0 = 0;
subject to v9_4_4:
	-z9_4_3 - y9_3_4 - x8_4_4 + x9_4_4 + y9_4_4 + z9_4_4 - 1.0 = 0;
subject to v10_4_4:
	-x9_4_4 - z10_4_3 - y10_3_4 + y10_4_4 + z10_4_4 + y11_4_4 + z11_4_4 - 1.0 = 0;
subject to v1_5_4:
	-z1_5_3 - y1_4_4 + x1_5_4 + y1_5_4 + z1_5_4 + y0_5_4 + z0_5_4 - 1.0 = 0;
subject to v2_5_4:
	-z2_5_3 - y2_4_4 - x1_5_4 + x2_5_4 + y2_5_4 + z2_5_4 - 1.0 = 0;
subject to v3_5_4:
	-z3_5_3 - y3_4_4 - x2_5_4 + x3_5_4 + y3_5_4 + z3_5_4 - 1.0 = 0;
subject to v4_5_4:
	-z4_5_3 - y4_4_4 - x3_5_4 + x4_5_4 + y4_5_4 + z4_5_4 - 1.0 = 0;
subject to v5_5_4:
	-z5_5_3 - y5_4_4 - x4_5_4 + x5_5_4 + y5_5_4 + z5_5_4 - 1.0 = 0;
subject to v6_5_4:
	-z6_5_3 - y6_4_4 - x5_5_4 + x6_5_4 + y6_5_4 + z6_5_4 - 1.0 = 0;
subject to v7_5_4:
	-z7_5_3 - y7_4_4 - x6_5_4 + x7_5_4 + y7_5_4 + z7_5_4 - 1.0 = 0;
subject to v8_5_4:
	-z8_5_3 - y8_4_4 - x7_5_4 + x8_5_4 + y8_5_4 + z8_5_4 - 1.0 = 0;
subject to v9_5_4:
	-z9_5_3 - y9_4_4 - x8_5_4 + x9_5_4 + y9_5_4 + z9_5_4 - 1.0 = 0;
subject to v10_5_4:
	-x9_5_4 - z10_5_3 - y10_4_4 + y10_5_4 + z10_5_4 + y11_5_4 + z11_5_4 - 1.0 = 0;
subject to v1_6_4:
	-z1_6_3 - y1_5_4 + x1_6_4 + y1_6_4 + z1_6_4 + y0_6_4 + z0_6_4 - 1.0 = 0;
subject to v2_6_4:
	-z2_6_3 - y2_5_4 - x1_6_4 + x2_6_4 + y2_6_4 + z2_6_4 - 1.0 = 0;
subject to v3_6_4:
	-z3_6_3 - y3_5_4 - x2_6_4 + x3_6_4 + y3_6_4 + z3_6_4 - 1.0 = 0;
subject to v4_6_4:
	-z4_6_3 - y4_5_4 - x3_6_4 + x4_6_4 + y4_6_4 + z4_6_4 - 1.0 = 0;
subject to v5_6_4:
	-z5_6_3 - y5_5_4 - x4_6_4 + x5_6_4 + y5_6_4 + z5_6_4 - 1.0 = 0;
subject to v6_6_4:
	-z6_6_3 - y6_5_4 - x5_6_4 + x6_6_4 + y6_6_4 + z6_6_4 - 1.0 = 0;
subject to v7_6_4:
	-z7_6_3 - y7_5_4 - x6_6_4 + x7_6_4 + y7_6_4 + z7_6_4 - 1.0 = 0;
subject to v8_6_4:
	-z8_6_3 - y8_5_4 - x7_6_4 + x8_6_4 + y8_6_4 + z8_6_4 - 1.0 = 0;
subject to v9_6_4:
	-z9_6_3 - y9_5_4 - x8_6_4 + x9_6_4 + y9_6_4 + z9_6_4 - 1.0 = 0;
subject to v10_6_4:
	-x9_6_4 - z10_6_3 - y10_5_4 + y10_6_4 + z10_6_4 + y11_6_4 + z11_6_4 - 1.0 = 0;
subject to v1_7_4:
	-z1_7_3 - y1_6_4 + x1_7_4 + y1_7_4 + z1_7_4 + y0_7_4 + z0_7_4 - 1.0 = 0;
subject to v2_7_4:
	-z2_7_3 - y2_6_4 - x1_7_4 + x2_7_4 + y2_7_4 + z2_7_4 - 1.0 = 0;
subject to v3_7_4:
	-z3_7_3 - y3_6_4 - x2_7_4 + x3_7_4 + y3_7_4 + z3_7_4 - 1.0 = 0;
subject to v4_7_4:
	-z4_7_3 - y4_6_4 - x3_7_4 + x4_7_4 + y4_7_4 + z4_7_4 - 1.0 = 0;
subject to v5_7_4:
	-z5_7_3 - y5_6_4 - x4_7_4 + x5_7_4 + y5_7_4 + z5_7_4 - 1.0 = 0;
subject to v6_7_4:
	-z6_7_3 - y6_6_4 - x5_7_4 + x6_7_4 + y6_7_4 + z6_7_4 - 1.0 = 0;
subject to v7_7_4:
	-z7_7_3 - y7_6_4 - x6_7_4 + x7_7_4 + y7_7_4 + z7_7_4 - 1.0 = 0;
subject to v8_7_4:
	-z8_7_3 - y8_6_4 - x7_7_4 + x8_7_4 + y8_7_4 + z8_7_4 - 1.0 = 0;
subject to v9_7_4:
	-z9_7_3 - y9_6_4 - x8_7_4 + x9_7_4 + y9_7_4 + z9_7_4 - 1.0 = 0;
subject to v10_7_4:
	-x9_7_4 - z10_7_3 - y10_6_4 + y10_7_4 + z10_7_4 + y11_7_4 + z11_7_4 - 1.0 = 0;
subject to v1_8_4:
	-z1_8_3 - y1_7_4 + x1_8_4 + y1_8_4 + z1_8_4 + y0_8_4 + z0_8_4 - 1.0 = 0;
subject to v2_8_4:
	-z2_8_3 - y2_7_4 - x1_8_4 + x2_8_4 + y2_8_4 + z2_8_4 - 1.0 = 0;
subject to v3_8_4:
	-z3_8_3 - y3_7_4 - x2_8_4 + x3_8_4 + y3_8_4 + z3_8_4 - 1.0 = 0;
subject to v4_8_4:
	-z4_8_3 - y4_7_4 - x3_8_4 + x4_8_4 + y4_8_4 + z4_8_4 - 1.0 = 0;
subject to v5_8_4:
	-z5_8_3 - y5_7_4 - x4_8_4 + x5_8_4 + y5_8_4 + z5_8_4 - 1.0 = 0;
subject to v6_8_4:
	-z6_8_3 - y6_7_4 - x5_8_4 + x6_8_4 + y6_8_4 + z6_8_4 - 1.0 = 0;
subject to v7_8_4:
	-z7_8_3 - y7_7_4 - x6_8_4 + x7_8_4 + y7_8_4 + z7_8_4 - 1.0 = 0;
subject to v8_8_4:
	-z8_8_3 - y8_7_4 - x7_8_4 + x8_8_4 + y8_8_4 + z8_8_4 - 1.0 = 0;
subject to v9_8_4:
	-z9_8_3 - y9_7_4 - x8_8_4 + x9_8_4 + y9_8_4 + z9_8_4 - 1.0 = 0;
subject to v10_8_4:
	-x9_8_4 - z10_8_3 - y10_7_4 + y10_8_4 + z10_8_4 + y11_8_4 + z11_8_4 - 1.0 = 0;
subject to v1_9_4:
	-z1_9_3 - y1_8_4 + x1_9_4 + y1_9_4 + z1_9_4 + y0_9_4 + z0_9_4 - 1.0 = 0;
subject to v2_9_4:
	-z2_9_3 - y2_8_4 - x1_9_4 + x2_9_4 + y2_9_4 + z2_9_4 - 1.0 = 0;
subject to v3_9_4:
	-z3_9_3 - y3_8_4 - x2_9_4 + x3_9_4 + y3_9_4 + z3_9_4 - 1.0 = 0;
subject to v4_9_4:
	-z4_9_3 - y4_8_4 - x3_9_4 + x4_9_4 + y4_9_4 + z4_9_4 - 1.0 = 0;
subject to v5_9_4:
	-z5_9_3 - y5_8_4 - x4_9_4 + x5_9_4 + y5_9_4 + z5_9_4 - 1.0 = 0;
subject to v6_9_4:
	-z6_9_3 - y6_8_4 - x5_9_4 + x6_9_4 + y6_9_4 + z6_9_4 - 1.0 = 0;
subject to v7_9_4:
	-z7_9_3 - y7_8_4 - x6_9_4 + x7_9_4 + y7_9_4 + z7_9_4 - 1.0 = 0;
subject to v8_9_4:
	-z8_9_3 - y8_8_4 - x7_9_4 + x8_9_4 + y8_9_4 + z8_9_4 - 1.0 = 0;
subject to v9_9_4:
	-z9_9_3 - y9_8_4 - x8_9_4 + x9_9_4 + y9_9_4 + z9_9_4 - 1.0 = 0;
subject to v10_9_4:
	-x9_9_4 - z10_9_3 - y10_8_4 + y10_9_4 + z10_9_4 + y11_9_4 + z11_9_4 - 1.0 = 0;
subject to v1_10_4:
	-y1_9_4 - z1_10_3 + x1_10_4 + z1_10_4 + y0_10_4 + z0_10_4 + x1_11_4 + z1_11_4 - 
	1.0 = 0;
subject to v2_10_4:
	-y2_9_4 - z2_10_3 - x1_10_4 + x2_10_4 + z2_10_4 + x2_11_4 + z2_11_4 - 1.0 = 0;
subject to v3_10_4:
	-y3_9_4 - z3_10_3 - x2_10_4 + x3_10_4 + z3_10_4 + x3_11_4 + z3_11_4 - 1.0 = 0;
subject to v4_10_4:
	-y4_9_4 - z4_10_3 - x3_10_4 + x4_10_4 + z4_10_4 + x4_11_4 + z4_11_4 - 1.0 = 0;
subject to v5_10_4:
	-y5_9_4 - z5_10_3 - x4_10_4 + x5_10_4 + z5_10_4 + x5_11_4 + z5_11_4 - 1.0 = 0;
subject to v6_10_4:
	-y6_9_4 - z6_10_3 - x5_10_4 + x6_10_4 + z6_10_4 + x6_11_4 + z6_11_4 - 1.0 = 0;
subject to v7_10_4:
	-y7_9_4 - z7_10_3 - x6_10_4 + x7_10_4 + z7_10_4 + x7_11_4 + z7_11_4 - 1.0 = 0;
subject to v8_10_4:
	-y8_9_4 - z8_10_3 - x7_10_4 + x8_10_4 + z8_10_4 + x8_11_4 + z8_11_4 - 1.0 = 0;
subject to v9_10_4:
	-y9_9_4 - z9_10_3 - x8_10_4 + x9_10_4 + z9_10_4 + x9_11_4 + z9_11_4 - 1.0 = 0;
subject to v10_10_4:
	-y10_9_4 - x9_10_4 + y11_10_4 + z11_10_4 + x10_11_4 + z10_11_4 - 1.0 = 0;
subject to v1_1_5:
	-z1_1_4 + x1_1_5 + y1_1_5 + z1_1_5 + y0_1_5 + z0_1_5 + x1_0_5 + z1_0_5 - 1.0 = 
	0;
subject to v2_1_5:
	-z2_1_4 - x1_1_5 + x2_1_5 + y2_1_5 + z2_1_5 + x2_0_5 + z2_0_5 - 1.0 = 0;
subject to v3_1_5:
	-z3_1_4 - x2_1_5 + x3_1_5 + y3_1_5 + z3_1_5 + x3_0_5 + z3_0_5 - 1.0 = 0;
subject to v4_1_5:
	-z4_1_4 - x3_1_5 + x4_1_5 + y4_1_5 + z4_1_5 + x4_0_5 + z4_0_5 - 1.0 = 0;
subject to v5_1_5:
	-z5_1_4 - x4_1_5 + x5_1_5 + y5_1_5 + z5_1_5 + x5_0_5 + z5_0_5 - 1.0 = 0;
subject to v6_1_5:
	-z6_1_4 - x5_1_5 + x6_1_5 + y6_1_5 + z6_1_5 + x6_0_5 + z6_0_5 - 1.0 = 0;
subject to v7_1_5:
	-z7_1_4 - x6_1_5 + x7_1_5 + y7_1_5 + z7_1_5 + x7_0_5 + z7_0_5 - 1.0 = 0;
subject to v8_1_5:
	-z8_1_4 - x7_1_5 + x8_1_5 + y8_1_5 + z8_1_5 + x8_0_5 + z8_0_5 - 1.0 = 0;
subject to v9_1_5:
	-z9_1_4 - x8_1_5 + x9_1_5 + y9_1_5 + z9_1_5 + x9_0_5 + z9_0_5 - 1.0 = 0;
subject to v10_1_5:
	-x9_1_5 - z10_1_4 + y10_1_5 + z10_1_5 + y11_1_5 + z11_1_5 + x10_0_5 + z10_0_5 - 
	1.0 = 0;
subject to v1_2_5:
	-z1_2_4 - y1_1_5 + x1_2_5 + y1_2_5 + z1_2_5 + y0_2_5 + z0_2_5 - 1.0 = 0;
subject to v2_2_5:
	-z2_2_4 - y2_1_5 - x1_2_5 + x2_2_5 + y2_2_5 + z2_2_5 - 1.0 = 0;
subject to v3_2_5:
	-z3_2_4 - y3_1_5 - x2_2_5 + x3_2_5 + y3_2_5 + z3_2_5 - 1.0 = 0;
subject to v4_2_5:
	-z4_2_4 - y4_1_5 - x3_2_5 + x4_2_5 + y4_2_5 + z4_2_5 - 1.0 = 0;
subject to v5_2_5:
	-z5_2_4 - y5_1_5 - x4_2_5 + x5_2_5 + y5_2_5 + z5_2_5 - 1.0 = 0;
subject to v6_2_5:
	-z6_2_4 - y6_1_5 - x5_2_5 + x6_2_5 + y6_2_5 + z6_2_5 - 1.0 = 0;
subject to v7_2_5:
	-z7_2_4 - y7_1_5 - x6_2_5 + x7_2_5 + y7_2_5 + z7_2_5 - 1.0 = 0;
subject to v8_2_5:
	-z8_2_4 - y8_1_5 - x7_2_5 + x8_2_5 + y8_2_5 + z8_2_5 - 1.0 = 0;
subject to v9_2_5:
	-z9_2_4 - y9_1_5 - x8_2_5 + x9_2_5 + y9_2_5 + z9_2_5 - 1.0 = 0;
subject to v10_2_5:
	-x9_2_5 - z10_2_4 - y10_1_5 + y10_2_5 + z10_2_5 + y11_2_5 + z11_2_5 - 1.0 = 0;
subject to v1_3_5:
	-z1_3_4 - y1_2_5 + x1_3_5 + y1_3_5 + z1_3_5 + y0_3_5 + z0_3_5 - 1.0 = 0;
subject to v2_3_5:
	-z2_3_4 - y2_2_5 - x1_3_5 + x2_3_5 + y2_3_5 + z2_3_5 - 1.0 = 0;
subject to v3_3_5:
	-z3_3_4 - y3_2_5 - x2_3_5 + x3_3_5 + y3_3_5 + z3_3_5 - 1.0 = 0;
subject to v4_3_5:
	-z4_3_4 - y4_2_5 - x3_3_5 + x4_3_5 + y4_3_5 + z4_3_5 - 1.0 = 0;
subject to v5_3_5:
	-z5_3_4 - y5_2_5 - x4_3_5 + x5_3_5 + y5_3_5 + z5_3_5 - 1.0 = 0;
subject to v6_3_5:
	-z6_3_4 - y6_2_5 - x5_3_5 + x6_3_5 + y6_3_5 + z6_3_5 - 1.0 = 0;
subject to v7_3_5:
	-z7_3_4 - y7_2_5 - x6_3_5 + x7_3_5 + y7_3_5 + z7_3_5 - 1.0 = 0;
subject to v8_3_5:
	-z8_3_4 - y8_2_5 - x7_3_5 + x8_3_5 + y8_3_5 + z8_3_5 - 1.0 = 0;
subject to v9_3_5:
	-z9_3_4 - y9_2_5 - x8_3_5 + x9_3_5 + y9_3_5 + z9_3_5 - 1.0 = 0;
subject to v10_3_5:
	-x9_3_5 - z10_3_4 - y10_2_5 + y10_3_5 + z10_3_5 + y11_3_5 + z11_3_5 - 1.0 = 0;
subject to v1_4_5:
	-z1_4_4 - y1_3_5 + x1_4_5 + y1_4_5 + z1_4_5 + y0_4_5 + z0_4_5 - 1.0 = 0;
subject to v2_4_5:
	-z2_4_4 - y2_3_5 - x1_4_5 + x2_4_5 + y2_4_5 + z2_4_5 - 1.0 = 0;
subject to v3_4_5:
	-z3_4_4 - y3_3_5 - x2_4_5 + x3_4_5 + y3_4_5 + z3_4_5 - 1.0 = 0;
subject to v4_4_5:
	-z4_4_4 - y4_3_5 - x3_4_5 + x4_4_5 + y4_4_5 + z4_4_5 - 1.0 = 0;
subject to v5_4_5:
	-z5_4_4 - y5_3_5 - x4_4_5 + x5_4_5 + y5_4_5 + z5_4_5 - 1.0 = 0;
subject to v6_4_5:
	-z6_4_4 - y6_3_5 - x5_4_5 + x6_4_5 + y6_4_5 + z6_4_5 - 1.0 = 0;
subject to v7_4_5:
	-z7_4_4 - y7_3_5 - x6_4_5 + x7_4_5 + y7_4_5 + z7_4_5 - 1.0 = 0;
subject to v8_4_5:
	-z8_4_4 - y8_3_5 - x7_4_5 + x8_4_5 + y8_4_5 + z8_4_5 - 1.0 = 0;
subject to v9_4_5:
	-z9_4_4 - y9_3_5 - x8_4_5 + x9_4_5 + y9_4_5 + z9_4_5 - 1.0 = 0;
subject to v10_4_5:
	-x9_4_5 - z10_4_4 - y10_3_5 + y10_4_5 + z10_4_5 + y11_4_5 + z11_4_5 - 1.0 = 0;
subject to v1_5_5:
	-z1_5_4 - y1_4_5 + x1_5_5 + y1_5_5 + z1_5_5 + y0_5_5 + z0_5_5 - 1.0 = 0;
subject to v2_5_5:
	-z2_5_4 - y2_4_5 - x1_5_5 + x2_5_5 + y2_5_5 + z2_5_5 - 1.0 = 0;
subject to v3_5_5:
	-z3_5_4 - y3_4_5 - x2_5_5 + x3_5_5 + y3_5_5 + z3_5_5 - 1.0 = 0;
subject to v4_5_5:
	-z4_5_4 - y4_4_5 - x3_5_5 + x4_5_5 + y4_5_5 + z4_5_5 - 1.0 = 0;
subject to v5_5_5:
	-z5_5_4 - y5_4_5 - x4_5_5 + x5_5_5 + y5_5_5 + z5_5_5 - 1.0 = 0;
subject to v6_5_5:
	-z6_5_4 - y6_4_5 - x5_5_5 + x6_5_5 + y6_5_5 + z6_5_5 - 1.0 = 0;
subject to v7_5_5:
	-z7_5_4 - y7_4_5 - x6_5_5 + x7_5_5 + y7_5_5 + z7_5_5 - 1.0 = 0;
subject to v8_5_5:
	-z8_5_4 - y8_4_5 - x7_5_5 + x8_5_5 + y8_5_5 + z8_5_5 - 1.0 = 0;
subject to v9_5_5:
	-z9_5_4 - y9_4_5 - x8_5_5 + x9_5_5 + y9_5_5 + z9_5_5 - 1.0 = 0;
subject to v10_5_5:
	-x9_5_5 - z10_5_4 - y10_4_5 + y10_5_5 + z10_5_5 + y11_5_5 + z11_5_5 - 1.0 = 0;
subject to v1_6_5:
	-z1_6_4 - y1_5_5 + x1_6_5 + y1_6_5 + z1_6_5 + y0_6_5 + z0_6_5 - 1.0 = 0;
subject to v2_6_5:
	-z2_6_4 - y2_5_5 - x1_6_5 + x2_6_5 + y2_6_5 + z2_6_5 - 1.0 = 0;
subject to v3_6_5:
	-z3_6_4 - y3_5_5 - x2_6_5 + x3_6_5 + y3_6_5 + z3_6_5 - 1.0 = 0;
subject to v4_6_5:
	-z4_6_4 - y4_5_5 - x3_6_5 + x4_6_5 + y4_6_5 + z4_6_5 - 1.0 = 0;
subject to v5_6_5:
	-z5_6_4 - y5_5_5 - x4_6_5 + x5_6_5 + y5_6_5 + z5_6_5 - 1.0 = 0;
subject to v6_6_5:
	-z6_6_4 - y6_5_5 - x5_6_5 + x6_6_5 + y6_6_5 + z6_6_5 - 1.0 = 0;
subject to v7_6_5:
	-z7_6_4 - y7_5_5 - x6_6_5 + x7_6_5 + y7_6_5 + z7_6_5 - 1.0 = 0;
subject to v8_6_5:
	-z8_6_4 - y8_5_5 - x7_6_5 + x8_6_5 + y8_6_5 + z8_6_5 - 1.0 = 0;
subject to v9_6_5:
	-z9_6_4 - y9_5_5 - x8_6_5 + x9_6_5 + y9_6_5 + z9_6_5 - 1.0 = 0;
subject to v10_6_5:
	-x9_6_5 - z10_6_4 - y10_5_5 + y10_6_5 + z10_6_5 + y11_6_5 + z11_6_5 - 1.0 = 0;
subject to v1_7_5:
	-z1_7_4 - y1_6_5 + x1_7_5 + y1_7_5 + z1_7_5 + y0_7_5 + z0_7_5 - 1.0 = 0;
subject to v2_7_5:
	-z2_7_4 - y2_6_5 - x1_7_5 + x2_7_5 + y2_7_5 + z2_7_5 - 1.0 = 0;
subject to v3_7_5:
	-z3_7_4 - y3_6_5 - x2_7_5 + x3_7_5 + y3_7_5 + z3_7_5 - 1.0 = 0;
subject to v4_7_5:
	-z4_7_4 - y4_6_5 - x3_7_5 + x4_7_5 + y4_7_5 + z4_7_5 - 1.0 = 0;
subject to v5_7_5:
	-z5_7_4 - y5_6_5 - x4_7_5 + x5_7_5 + y5_7_5 + z5_7_5 - 1.0 = 0;
subject to v6_7_5:
	-z6_7_4 - y6_6_5 - x5_7_5 + x6_7_5 + y6_7_5 + z6_7_5 - 1.0 = 0;
subject to v7_7_5:
	-z7_7_4 - y7_6_5 - x6_7_5 + x7_7_5 + y7_7_5 + z7_7_5 - 1.0 = 0;
subject to v8_7_5:
	-z8_7_4 - y8_6_5 - x7_7_5 + x8_7_5 + y8_7_5 + z8_7_5 - 1.0 = 0;
subject to v9_7_5:
	-z9_7_4 - y9_6_5 - x8_7_5 + x9_7_5 + y9_7_5 + z9_7_5 - 1.0 = 0;
subject to v10_7_5:
	-x9_7_5 - z10_7_4 - y10_6_5 + y10_7_5 + z10_7_5 + y11_7_5 + z11_7_5 - 1.0 = 0;
subject to v1_8_5:
	-z1_8_4 - y1_7_5 + x1_8_5 + y1_8_5 + z1_8_5 + y0_8_5 + z0_8_5 - 1.0 = 0;
subject to v2_8_5:
	-z2_8_4 - y2_7_5 - x1_8_5 + x2_8_5 + y2_8_5 + z2_8_5 - 1.0 = 0;
subject to v3_8_5:
	-z3_8_4 - y3_7_5 - x2_8_5 + x3_8_5 + y3_8_5 + z3_8_5 - 1.0 = 0;
subject to v4_8_5:
	-z4_8_4 - y4_7_5 - x3_8_5 + x4_8_5 + y4_8_5 + z4_8_5 - 1.0 = 0;
subject to v5_8_5:
	-z5_8_4 - y5_7_5 - x4_8_5 + x5_8_5 + y5_8_5 + z5_8_5 - 1.0 = 0;
subject to v6_8_5:
	-z6_8_4 - y6_7_5 - x5_8_5 + x6_8_5 + y6_8_5 + z6_8_5 - 1.0 = 0;
subject to v7_8_5:
	-z7_8_4 - y7_7_5 - x6_8_5 + x7_8_5 + y7_8_5 + z7_8_5 - 1.0 = 0;
subject to v8_8_5:
	-z8_8_4 - y8_7_5 - x7_8_5 + x8_8_5 + y8_8_5 + z8_8_5 - 1.0 = 0;
subject to v9_8_5:
	-z9_8_4 - y9_7_5 - x8_8_5 + x9_8_5 + y9_8_5 + z9_8_5 - 1.0 = 0;
subject to v10_8_5:
	-x9_8_5 - z10_8_4 - y10_7_5 + y10_8_5 + z10_8_5 + y11_8_5 + z11_8_5 - 1.0 = 0;
subject to v1_9_5:
	-z1_9_4 - y1_8_5 + x1_9_5 + y1_9_5 + z1_9_5 + y0_9_5 + z0_9_5 - 1.0 = 0;
subject to v2_9_5:
	-z2_9_4 - y2_8_5 - x1_9_5 + x2_9_5 + y2_9_5 + z2_9_5 - 1.0 = 0;
subject to v3_9_5:
	-z3_9_4 - y3_8_5 - x2_9_5 + x3_9_5 + y3_9_5 + z3_9_5 - 1.0 = 0;
subject to v4_9_5:
	-z4_9_4 - y4_8_5 - x3_9_5 + x4_9_5 + y4_9_5 + z4_9_5 - 1.0 = 0;
subject to v5_9_5:
	-z5_9_4 - y5_8_5 - x4_9_5 + x5_9_5 + y5_9_5 + z5_9_5 - 1.0 = 0;
subject to v6_9_5:
	-z6_9_4 - y6_8_5 - x5_9_5 + x6_9_5 + y6_9_5 + z6_9_5 - 1.0 = 0;
subject to v7_9_5:
	-z7_9_4 - y7_8_5 - x6_9_5 + x7_9_5 + y7_9_5 + z7_9_5 - 1.0 = 0;
subject to v8_9_5:
	-z8_9_4 - y8_8_5 - x7_9_5 + x8_9_5 + y8_9_5 + z8_9_5 - 1.0 = 0;
subject to v9_9_5:
	-z9_9_4 - y9_8_5 - x8_9_5 + x9_9_5 + y9_9_5 + z9_9_5 - 1.0 = 0;
subject to v10_9_5:
	-x9_9_5 - z10_9_4 - y10_8_5 + y10_9_5 + z10_9_5 + y11_9_5 + z11_9_5 - 1.0 = 0;
subject to v1_10_5:
	-y1_9_5 - z1_10_4 + x1_10_5 + z1_10_5 + y0_10_5 + z0_10_5 + x1_11_5 + z1_11_5 - 
	1.0 = 0;
subject to v2_10_5:
	-y2_9_5 - z2_10_4 - x1_10_5 + x2_10_5 + z2_10_5 + x2_11_5 + z2_11_5 - 1.0 = 0;
subject to v3_10_5:
	-y3_9_5 - z3_10_4 - x2_10_5 + x3_10_5 + z3_10_5 + x3_11_5 + z3_11_5 - 1.0 = 0;
subject to v4_10_5:
	-y4_9_5 - z4_10_4 - x3_10_5 + x4_10_5 + z4_10_5 + x4_11_5 + z4_11_5 - 1.0 = 0;
subject to v5_10_5:
	-y5_9_5 - z5_10_4 - x4_10_5 + x5_10_5 + z5_10_5 + x5_11_5 + z5_11_5 - 1.0 = 0;
subject to v6_10_5:
	-y6_9_5 - z6_10_4 - x5_10_5 + x6_10_5 + z6_10_5 + x6_11_5 + z6_11_5 - 1.0 = 0;
subject to v7_10_5:
	-y7_9_5 - z7_10_4 - x6_10_5 + x7_10_5 + z7_10_5 + x7_11_5 + z7_11_5 - 1.0 = 0;
subject to v8_10_5:
	-y8_9_5 - z8_10_4 - x7_10_5 + x8_10_5 + z8_10_5 + x8_11_5 + z8_11_5 - 1.0 = 0;
subject to v9_10_5:
	-y9_9_5 - z9_10_4 - x8_10_5 + x9_10_5 + z9_10_5 + x9_11_5 + z9_11_5 - 1.0 = 0;
subject to v10_10_5:
	-y10_9_5 - x9_10_5 + y11_10_5 + z11_10_5 + x10_11_5 + z10_11_5 - 1.0 = 0;
subject to v1_1_6:
	-z1_1_5 + x1_1_6 + y1_1_6 + z1_1_6 + y0_1_6 + z0_1_6 + x1_0_6 + z1_0_6 - 1.0 = 
	0;
subject to v2_1_6:
	-z2_1_5 - x1_1_6 + x2_1_6 + y2_1_6 + z2_1_6 + x2_0_6 + z2_0_6 - 1.0 = 0;
subject to v3_1_6:
	-z3_1_5 - x2_1_6 + x3_1_6 + y3_1_6 + z3_1_6 + x3_0_6 + z3_0_6 - 1.0 = 0;
subject to v4_1_6:
	-z4_1_5 - x3_1_6 + x4_1_6 + y4_1_6 + z4_1_6 + x4_0_6 + z4_0_6 - 1.0 = 0;
subject to v5_1_6:
	-z5_1_5 - x4_1_6 + x5_1_6 + y5_1_6 + z5_1_6 + x5_0_6 + z5_0_6 - 1.0 = 0;
subject to v6_1_6:
	-z6_1_5 - x5_1_6 + x6_1_6 + y6_1_6 + z6_1_6 + x6_0_6 + z6_0_6 - 1.0 = 0;
subject to v7_1_6:
	-z7_1_5 - x6_1_6 + x7_1_6 + y7_1_6 + z7_1_6 + x7_0_6 + z7_0_6 - 1.0 = 0;
subject to v8_1_6:
	-z8_1_5 - x7_1_6 + x8_1_6 + y8_1_6 + z8_1_6 + x8_0_6 + z8_0_6 - 1.0 = 0;
subject to v9_1_6:
	-z9_1_5 - x8_1_6 + x9_1_6 + y9_1_6 + z9_1_6 + x9_0_6 + z9_0_6 - 1.0 = 0;
subject to v10_1_6:
	-x9_1_6 - z10_1_5 + y10_1_6 + z10_1_6 + y11_1_6 + z11_1_6 + x10_0_6 + z10_0_6 - 
	1.0 = 0;
subject to v1_2_6:
	-z1_2_5 - y1_1_6 + x1_2_6 + y1_2_6 + z1_2_6 + y0_2_6 + z0_2_6 - 1.0 = 0;
subject to v2_2_6:
	-z2_2_5 - y2_1_6 - x1_2_6 + x2_2_6 + y2_2_6 + z2_2_6 - 1.0 = 0;
subject to v3_2_6:
	-z3_2_5 - y3_1_6 - x2_2_6 + x3_2_6 + y3_2_6 + z3_2_6 - 1.0 = 0;
subject to v4_2_6:
	-z4_2_5 - y4_1_6 - x3_2_6 + x4_2_6 + y4_2_6 + z4_2_6 - 1.0 = 0;
subject to v5_2_6:
	-z5_2_5 - y5_1_6 - x4_2_6 + x5_2_6 + y5_2_6 + z5_2_6 - 1.0 = 0;
subject to v6_2_6:
	-z6_2_5 - y6_1_6 - x5_2_6 + x6_2_6 + y6_2_6 + z6_2_6 - 1.0 = 0;
subject to v7_2_6:
	-z7_2_5 - y7_1_6 - x6_2_6 + x7_2_6 + y7_2_6 + z7_2_6 - 1.0 = 0;
subject to v8_2_6:
	-z8_2_5 - y8_1_6 - x7_2_6 + x8_2_6 + y8_2_6 + z8_2_6 - 1.0 = 0;
subject to v9_2_6:
	-z9_2_5 - y9_1_6 - x8_2_6 + x9_2_6 + y9_2_6 + z9_2_6 - 1.0 = 0;
subject to v10_2_6:
	-x9_2_6 - z10_2_5 - y10_1_6 + y10_2_6 + z10_2_6 + y11_2_6 + z11_2_6 - 1.0 = 0;
subject to v1_3_6:
	-z1_3_5 - y1_2_6 + x1_3_6 + y1_3_6 + z1_3_6 + y0_3_6 + z0_3_6 - 1.0 = 0;
subject to v2_3_6:
	-z2_3_5 - y2_2_6 - x1_3_6 + x2_3_6 + y2_3_6 + z2_3_6 - 1.0 = 0;
subject to v3_3_6:
	-z3_3_5 - y3_2_6 - x2_3_6 + x3_3_6 + y3_3_6 + z3_3_6 - 1.0 = 0;
subject to v4_3_6:
	-z4_3_5 - y4_2_6 - x3_3_6 + x4_3_6 + y4_3_6 + z4_3_6 - 1.0 = 0;
subject to v5_3_6:
	-z5_3_5 - y5_2_6 - x4_3_6 + x5_3_6 + y5_3_6 + z5_3_6 - 1.0 = 0;
subject to v6_3_6:
	-z6_3_5 - y6_2_6 - x5_3_6 + x6_3_6 + y6_3_6 + z6_3_6 - 1.0 = 0;
subject to v7_3_6:
	-z7_3_5 - y7_2_6 - x6_3_6 + x7_3_6 + y7_3_6 + z7_3_6 - 1.0 = 0;
subject to v8_3_6:
	-z8_3_5 - y8_2_6 - x7_3_6 + x8_3_6 + y8_3_6 + z8_3_6 - 1.0 = 0;
subject to v9_3_6:
	-z9_3_5 - y9_2_6 - x8_3_6 + x9_3_6 + y9_3_6 + z9_3_6 - 1.0 = 0;
subject to v10_3_6:
	-x9_3_6 - z10_3_5 - y10_2_6 + y10_3_6 + z10_3_6 + y11_3_6 + z11_3_6 - 1.0 = 0;
subject to v1_4_6:
	-z1_4_5 - y1_3_6 + x1_4_6 + y1_4_6 + z1_4_6 + y0_4_6 + z0_4_6 - 1.0 = 0;
subject to v2_4_6:
	-z2_4_5 - y2_3_6 - x1_4_6 + x2_4_6 + y2_4_6 + z2_4_6 - 1.0 = 0;
subject to v3_4_6:
	-z3_4_5 - y3_3_6 - x2_4_6 + x3_4_6 + y3_4_6 + z3_4_6 - 1.0 = 0;
subject to v4_4_6:
	-z4_4_5 - y4_3_6 - x3_4_6 + x4_4_6 + y4_4_6 + z4_4_6 - 1.0 = 0;
subject to v5_4_6:
	-z5_4_5 - y5_3_6 - x4_4_6 + x5_4_6 + y5_4_6 + z5_4_6 - 1.0 = 0;
subject to v6_4_6:
	-z6_4_5 - y6_3_6 - x5_4_6 + x6_4_6 + y6_4_6 + z6_4_6 - 1.0 = 0;
subject to v7_4_6:
	-z7_4_5 - y7_3_6 - x6_4_6 + x7_4_6 + y7_4_6 + z7_4_6 - 1.0 = 0;
subject to v8_4_6:
	-z8_4_5 - y8_3_6 - x7_4_6 + x8_4_6 + y8_4_6 + z8_4_6 - 1.0 = 0;
subject to v9_4_6:
	-z9_4_5 - y9_3_6 - x8_4_6 + x9_4_6 + y9_4_6 + z9_4_6 - 1.0 = 0;
subject to v10_4_6:
	-x9_4_6 - z10_4_5 - y10_3_6 + y10_4_6 + z10_4_6 + y11_4_6 + z11_4_6 - 1.0 = 0;
subject to v1_5_6:
	-z1_5_5 - y1_4_6 + x1_5_6 + y1_5_6 + z1_5_6 + y0_5_6 + z0_5_6 - 1.0 = 0;
subject to v2_5_6:
	-z2_5_5 - y2_4_6 - x1_5_6 + x2_5_6 + y2_5_6 + z2_5_6 - 1.0 = 0;
subject to v3_5_6:
	-z3_5_5 - y3_4_6 - x2_5_6 + x3_5_6 + y3_5_6 + z3_5_6 - 1.0 = 0;
subject to v4_5_6:
	-z4_5_5 - y4_4_6 - x3_5_6 + x4_5_6 + y4_5_6 + z4_5_6 - 1.0 = 0;
subject to v5_5_6:
	-z5_5_5 - y5_4_6 - x4_5_6 + x5_5_6 + y5_5_6 + z5_5_6 - 1.0 = 0;
subject to v6_5_6:
	-z6_5_5 - y6_4_6 - x5_5_6 + x6_5_6 + y6_5_6 + z6_5_6 - 1.0 = 0;
subject to v7_5_6:
	-z7_5_5 - y7_4_6 - x6_5_6 + x7_5_6 + y7_5_6 + z7_5_6 - 1.0 = 0;
subject to v8_5_6:
	-z8_5_5 - y8_4_6 - x7_5_6 + x8_5_6 + y8_5_6 + z8_5_6 - 1.0 = 0;
subject to v9_5_6:
	-z9_5_5 - y9_4_6 - x8_5_6 + x9_5_6 + y9_5_6 + z9_5_6 - 1.0 = 0;
subject to v10_5_6:
	-x9_5_6 - z10_5_5 - y10_4_6 + y10_5_6 + z10_5_6 + y11_5_6 + z11_5_6 - 1.0 = 0;
subject to v1_6_6:
	-z1_6_5 - y1_5_6 + x1_6_6 + y1_6_6 + z1_6_6 + y0_6_6 + z0_6_6 - 1.0 = 0;
subject to v2_6_6:
	-z2_6_5 - y2_5_6 - x1_6_6 + x2_6_6 + y2_6_6 + z2_6_6 - 1.0 = 0;
subject to v3_6_6:
	-z3_6_5 - y3_5_6 - x2_6_6 + x3_6_6 + y3_6_6 + z3_6_6 - 1.0 = 0;
subject to v4_6_6:
	-z4_6_5 - y4_5_6 - x3_6_6 + x4_6_6 + y4_6_6 + z4_6_6 - 1.0 = 0;
subject to v5_6_6:
	-z5_6_5 - y5_5_6 - x4_6_6 + x5_6_6 + y5_6_6 + z5_6_6 - 1.0 = 0;
subject to v6_6_6:
	-z6_6_5 - y6_5_6 - x5_6_6 + x6_6_6 + y6_6_6 + z6_6_6 - 1.0 = 0;
subject to v7_6_6:
	-z7_6_5 - y7_5_6 - x6_6_6 + x7_6_6 + y7_6_6 + z7_6_6 - 1.0 = 0;
subject to v8_6_6:
	-z8_6_5 - y8_5_6 - x7_6_6 + x8_6_6 + y8_6_6 + z8_6_6 - 1.0 = 0;
subject to v9_6_6:
	-z9_6_5 - y9_5_6 - x8_6_6 + x9_6_6 + y9_6_6 + z9_6_6 - 1.0 = 0;
subject to v10_6_6:
	-x9_6_6 - z10_6_5 - y10_5_6 + y10_6_6 + z10_6_6 + y11_6_6 + z11_6_6 - 1.0 = 0;
subject to v1_7_6:
	-z1_7_5 - y1_6_6 + x1_7_6 + y1_7_6 + z1_7_6 + y0_7_6 + z0_7_6 - 1.0 = 0;
subject to v2_7_6:
	-z2_7_5 - y2_6_6 - x1_7_6 + x2_7_6 + y2_7_6 + z2_7_6 - 1.0 = 0;
subject to v3_7_6:
	-z3_7_5 - y3_6_6 - x2_7_6 + x3_7_6 + y3_7_6 + z3_7_6 - 1.0 = 0;
subject to v4_7_6:
	-z4_7_5 - y4_6_6 - x3_7_6 + x4_7_6 + y4_7_6 + z4_7_6 - 1.0 = 0;
subject to v5_7_6:
	-z5_7_5 - y5_6_6 - x4_7_6 + x5_7_6 + y5_7_6 + z5_7_6 - 1.0 = 0;
subject to v6_7_6:
	-z6_7_5 - y6_6_6 - x5_7_6 + x6_7_6 + y6_7_6 + z6_7_6 - 1.0 = 0;
subject to v7_7_6:
	-z7_7_5 - y7_6_6 - x6_7_6 + x7_7_6 + y7_7_6 + z7_7_6 - 1.0 = 0;
subject to v8_7_6:
	-z8_7_5 - y8_6_6 - x7_7_6 + x8_7_6 + y8_7_6 + z8_7_6 - 1.0 = 0;
subject to v9_7_6:
	-z9_7_5 - y9_6_6 - x8_7_6 + x9_7_6 + y9_7_6 + z9_7_6 - 1.0 = 0;
subject to v10_7_6:
	-x9_7_6 - z10_7_5 - y10_6_6 + y10_7_6 + z10_7_6 + y11_7_6 + z11_7_6 - 1.0 = 0;
subject to v1_8_6:
	-z1_8_5 - y1_7_6 + x1_8_6 + y1_8_6 + z1_8_6 + y0_8_6 + z0_8_6 - 1.0 = 0;
subject to v2_8_6:
	-z2_8_5 - y2_7_6 - x1_8_6 + x2_8_6 + y2_8_6 + z2_8_6 - 1.0 = 0;
subject to v3_8_6:
	-z3_8_5 - y3_7_6 - x2_8_6 + x3_8_6 + y3_8_6 + z3_8_6 - 1.0 = 0;
subject to v4_8_6:
	-z4_8_5 - y4_7_6 - x3_8_6 + x4_8_6 + y4_8_6 + z4_8_6 - 1.0 = 0;
subject to v5_8_6:
	-z5_8_5 - y5_7_6 - x4_8_6 + x5_8_6 + y5_8_6 + z5_8_6 - 1.0 = 0;
subject to v6_8_6:
	-z6_8_5 - y6_7_6 - x5_8_6 + x6_8_6 + y6_8_6 + z6_8_6 - 1.0 = 0;
subject to v7_8_6:
	-z7_8_5 - y7_7_6 - x6_8_6 + x7_8_6 + y7_8_6 + z7_8_6 - 1.0 = 0;
subject to v8_8_6:
	-z8_8_5 - y8_7_6 - x7_8_6 + x8_8_6 + y8_8_6 + z8_8_6 - 1.0 = 0;
subject to v9_8_6:
	-z9_8_5 - y9_7_6 - x8_8_6 + x9_8_6 + y9_8_6 + z9_8_6 - 1.0 = 0;
subject to v10_8_6:
	-x9_8_6 - z10_8_5 - y10_7_6 + y10_8_6 + z10_8_6 + y11_8_6 + z11_8_6 - 1.0 = 0;
subject to v1_9_6:
	-z1_9_5 - y1_8_6 + x1_9_6 + y1_9_6 + z1_9_6 + y0_9_6 + z0_9_6 - 1.0 = 0;
subject to v2_9_6:
	-z2_9_5 - y2_8_6 - x1_9_6 + x2_9_6 + y2_9_6 + z2_9_6 - 1.0 = 0;
subject to v3_9_6:
	-z3_9_5 - y3_8_6 - x2_9_6 + x3_9_6 + y3_9_6 + z3_9_6 - 1.0 = 0;
subject to v4_9_6:
	-z4_9_5 - y4_8_6 - x3_9_6 + x4_9_6 + y4_9_6 + z4_9_6 - 1.0 = 0;
subject to v5_9_6:
	-z5_9_5 - y5_8_6 - x4_9_6 + x5_9_6 + y5_9_6 + z5_9_6 - 1.0 = 0;
subject to v6_9_6:
	-z6_9_5 - y6_8_6 - x5_9_6 + x6_9_6 + y6_9_6 + z6_9_6 - 1.0 = 0;
subject to v7_9_6:
	-z7_9_5 - y7_8_6 - x6_9_6 + x7_9_6 + y7_9_6 + z7_9_6 - 1.0 = 0;
subject to v8_9_6:
	-z8_9_5 - y8_8_6 - x7_9_6 + x8_9_6 + y8_9_6 + z8_9_6 - 1.0 = 0;
subject to v9_9_6:
	-z9_9_5 - y9_8_6 - x8_9_6 + x9_9_6 + y9_9_6 + z9_9_6 - 1.0 = 0;
subject to v10_9_6:
	-x9_9_6 - z10_9_5 - y10_8_6 + y10_9_6 + z10_9_6 + y11_9_6 + z11_9_6 - 1.0 = 0;
subject to v1_10_6:
	-y1_9_6 - z1_10_5 + x1_10_6 + z1_10_6 + y0_10_6 + z0_10_6 + x1_11_6 + z1_11_6 - 
	1.0 = 0;
subject to v2_10_6:
	-y2_9_6 - z2_10_5 - x1_10_6 + x2_10_6 + z2_10_6 + x2_11_6 + z2_11_6 - 1.0 = 0;
subject to v3_10_6:
	-y3_9_6 - z3_10_5 - x2_10_6 + x3_10_6 + z3_10_6 + x3_11_6 + z3_11_6 - 1.0 = 0;
subject to v4_10_6:
	-y4_9_6 - z4_10_5 - x3_10_6 + x4_10_6 + z4_10_6 + x4_11_6 + z4_11_6 - 1.0 = 0;
subject to v5_10_6:
	-y5_9_6 - z5_10_5 - x4_10_6 + x5_10_6 + z5_10_6 + x5_11_6 + z5_11_6 - 1.0 = 0;
subject to v6_10_6:
	-y6_9_6 - z6_10_5 - x5_10_6 + x6_10_6 + z6_10_6 + x6_11_6 + z6_11_6 - 1.0 = 0;
subject to v7_10_6:
	-y7_9_6 - z7_10_5 - x6_10_6 + x7_10_6 + z7_10_6 + x7_11_6 + z7_11_6 - 1.0 = 0;
subject to v8_10_6:
	-y8_9_6 - z8_10_5 - x7_10_6 + x8_10_6 + z8_10_6 + x8_11_6 + z8_11_6 - 1.0 = 0;
subject to v9_10_6:
	-y9_9_6 - z9_10_5 - x8_10_6 + x9_10_6 + z9_10_6 + x9_11_6 + z9_11_6 - 1.0 = 0;
subject to v10_10_6:
	-y10_9_6 - x9_10_6 + y11_10_6 + z11_10_6 + x10_11_6 + z10_11_6 - 1.0 = 0;
subject to v1_1_7:
	-z1_1_6 + x1_1_7 + y1_1_7 + z1_1_7 + y0_1_7 + z0_1_7 + x1_0_7 + z1_0_7 - 1.0 = 
	0;
subject to v2_1_7:
	-z2_1_6 - x1_1_7 + x2_1_7 + y2_1_7 + z2_1_7 + x2_0_7 + z2_0_7 - 1.0 = 0;
subject to v3_1_7:
	-z3_1_6 - x2_1_7 + x3_1_7 + y3_1_7 + z3_1_7 + x3_0_7 + z3_0_7 - 1.0 = 0;
subject to v4_1_7:
	-z4_1_6 - x3_1_7 + x4_1_7 + y4_1_7 + z4_1_7 + x4_0_7 + z4_0_7 - 1.0 = 0;
subject to v5_1_7:
	-z5_1_6 - x4_1_7 + x5_1_7 + y5_1_7 + z5_1_7 + x5_0_7 + z5_0_7 - 1.0 = 0;
subject to v6_1_7:
	-z6_1_6 - x5_1_7 + x6_1_7 + y6_1_7 + z6_1_7 + x6_0_7 + z6_0_7 - 1.0 = 0;
subject to v7_1_7:
	-z7_1_6 - x6_1_7 + x7_1_7 + y7_1_7 + z7_1_7 + x7_0_7 + z7_0_7 - 1.0 = 0;
subject to v8_1_7:
	-z8_1_6 - x7_1_7 + x8_1_7 + y8_1_7 + z8_1_7 + x8_0_7 + z8_0_7 - 1.0 = 0;
subject to v9_1_7:
	-z9_1_6 - x8_1_7 + x9_1_7 + y9_1_7 + z9_1_7 + x9_0_7 + z9_0_7 - 1.0 = 0;
subject to v10_1_7:
	-x9_1_7 - z10_1_6 + y10_1_7 + z10_1_7 + y11_1_7 + z11_1_7 + x10_0_7 + z10_0_7 - 
	1.0 = 0;
subject to v1_2_7:
	-z1_2_6 - y1_1_7 + x1_2_7 + y1_2_7 + z1_2_7 + y0_2_7 + z0_2_7 - 1.0 = 0;
subject to v2_2_7:
	-z2_2_6 - y2_1_7 - x1_2_7 + x2_2_7 + y2_2_7 + z2_2_7 - 1.0 = 0;
subject to v3_2_7:
	-z3_2_6 - y3_1_7 - x2_2_7 + x3_2_7 + y3_2_7 + z3_2_7 - 1.0 = 0;
subject to v4_2_7:
	-z4_2_6 - y4_1_7 - x3_2_7 + x4_2_7 + y4_2_7 + z4_2_7 - 1.0 = 0;
subject to v5_2_7:
	-z5_2_6 - y5_1_7 - x4_2_7 + x5_2_7 + y5_2_7 + z5_2_7 - 1.0 = 0;
subject to v6_2_7:
	-z6_2_6 - y6_1_7 - x5_2_7 + x6_2_7 + y6_2_7 + z6_2_7 - 1.0 = 0;
subject to v7_2_7:
	-z7_2_6 - y7_1_7 - x6_2_7 + x7_2_7 + y7_2_7 + z7_2_7 - 1.0 = 0;
subject to v8_2_7:
	-z8_2_6 - y8_1_7 - x7_2_7 + x8_2_7 + y8_2_7 + z8_2_7 - 1.0 = 0;
subject to v9_2_7:
	-z9_2_6 - y9_1_7 - x8_2_7 + x9_2_7 + y9_2_7 + z9_2_7 - 1.0 = 0;
subject to v10_2_7:
	-x9_2_7 - z10_2_6 - y10_1_7 + y10_2_7 + z10_2_7 + y11_2_7 + z11_2_7 - 1.0 = 0;
subject to v1_3_7:
	-z1_3_6 - y1_2_7 + x1_3_7 + y1_3_7 + z1_3_7 + y0_3_7 + z0_3_7 - 1.0 = 0;
subject to v2_3_7:
	-z2_3_6 - y2_2_7 - x1_3_7 + x2_3_7 + y2_3_7 + z2_3_7 - 1.0 = 0;
subject to v3_3_7:
	-z3_3_6 - y3_2_7 - x2_3_7 + x3_3_7 + y3_3_7 + z3_3_7 - 1.0 = 0;
subject to v4_3_7:
	-z4_3_6 - y4_2_7 - x3_3_7 + x4_3_7 + y4_3_7 + z4_3_7 - 1.0 = 0;
subject to v5_3_7:
	-z5_3_6 - y5_2_7 - x4_3_7 + x5_3_7 + y5_3_7 + z5_3_7 - 1.0 = 0;
subject to v6_3_7:
	-z6_3_6 - y6_2_7 - x5_3_7 + x6_3_7 + y6_3_7 + z6_3_7 - 1.0 = 0;
subject to v7_3_7:
	-z7_3_6 - y7_2_7 - x6_3_7 + x7_3_7 + y7_3_7 + z7_3_7 - 1.0 = 0;
subject to v8_3_7:
	-z8_3_6 - y8_2_7 - x7_3_7 + x8_3_7 + y8_3_7 + z8_3_7 - 1.0 = 0;
subject to v9_3_7:
	-z9_3_6 - y9_2_7 - x8_3_7 + x9_3_7 + y9_3_7 + z9_3_7 - 1.0 = 0;
subject to v10_3_7:
	-x9_3_7 - z10_3_6 - y10_2_7 + y10_3_7 + z10_3_7 + y11_3_7 + z11_3_7 - 1.0 = 0;
subject to v1_4_7:
	-z1_4_6 - y1_3_7 + x1_4_7 + y1_4_7 + z1_4_7 + y0_4_7 + z0_4_7 - 1.0 = 0;
subject to v2_4_7:
	-z2_4_6 - y2_3_7 - x1_4_7 + x2_4_7 + y2_4_7 + z2_4_7 - 1.0 = 0;
subject to v3_4_7:
	-z3_4_6 - y3_3_7 - x2_4_7 + x3_4_7 + y3_4_7 + z3_4_7 - 1.0 = 0;
subject to v4_4_7:
	-z4_4_6 - y4_3_7 - x3_4_7 + x4_4_7 + y4_4_7 + z4_4_7 - 1.0 = 0;
subject to v5_4_7:
	-z5_4_6 - y5_3_7 - x4_4_7 + x5_4_7 + y5_4_7 + z5_4_7 - 1.0 = 0;
subject to v6_4_7:
	-z6_4_6 - y6_3_7 - x5_4_7 + x6_4_7 + y6_4_7 + z6_4_7 - 1.0 = 0;
subject to v7_4_7:
	-z7_4_6 - y7_3_7 - x6_4_7 + x7_4_7 + y7_4_7 + z7_4_7 - 1.0 = 0;
subject to v8_4_7:
	-z8_4_6 - y8_3_7 - x7_4_7 + x8_4_7 + y8_4_7 + z8_4_7 - 1.0 = 0;
subject to v9_4_7:
	-z9_4_6 - y9_3_7 - x8_4_7 + x9_4_7 + y9_4_7 + z9_4_7 - 1.0 = 0;
subject to v10_4_7:
	-x9_4_7 - z10_4_6 - y10_3_7 + y10_4_7 + z10_4_7 + y11_4_7 + z11_4_7 - 1.0 = 0;
subject to v1_5_7:
	-z1_5_6 - y1_4_7 + x1_5_7 + y1_5_7 + z1_5_7 + y0_5_7 + z0_5_7 - 1.0 = 0;
subject to v2_5_7:
	-z2_5_6 - y2_4_7 - x1_5_7 + x2_5_7 + y2_5_7 + z2_5_7 - 1.0 = 0;
subject to v3_5_7:
	-z3_5_6 - y3_4_7 - x2_5_7 + x3_5_7 + y3_5_7 + z3_5_7 - 1.0 = 0;
subject to v4_5_7:
	-z4_5_6 - y4_4_7 - x3_5_7 + x4_5_7 + y4_5_7 + z4_5_7 - 1.0 = 0;
subject to v5_5_7:
	-z5_5_6 - y5_4_7 - x4_5_7 + x5_5_7 + y5_5_7 + z5_5_7 - 1.0 = 0;
subject to v6_5_7:
	-z6_5_6 - y6_4_7 - x5_5_7 + x6_5_7 + y6_5_7 + z6_5_7 - 1.0 = 0;
subject to v7_5_7:
	-z7_5_6 - y7_4_7 - x6_5_7 + x7_5_7 + y7_5_7 + z7_5_7 - 1.0 = 0;
subject to v8_5_7:
	-z8_5_6 - y8_4_7 - x7_5_7 + x8_5_7 + y8_5_7 + z8_5_7 - 1.0 = 0;
subject to v9_5_7:
	-z9_5_6 - y9_4_7 - x8_5_7 + x9_5_7 + y9_5_7 + z9_5_7 - 1.0 = 0;
subject to v10_5_7:
	-x9_5_7 - z10_5_6 - y10_4_7 + y10_5_7 + z10_5_7 + y11_5_7 + z11_5_7 - 1.0 = 0;
subject to v1_6_7:
	-z1_6_6 - y1_5_7 + x1_6_7 + y1_6_7 + z1_6_7 + y0_6_7 + z0_6_7 - 1.0 = 0;
subject to v2_6_7:
	-z2_6_6 - y2_5_7 - x1_6_7 + x2_6_7 + y2_6_7 + z2_6_7 - 1.0 = 0;
subject to v3_6_7:
	-z3_6_6 - y3_5_7 - x2_6_7 + x3_6_7 + y3_6_7 + z3_6_7 - 1.0 = 0;
subject to v4_6_7:
	-z4_6_6 - y4_5_7 - x3_6_7 + x4_6_7 + y4_6_7 + z4_6_7 - 1.0 = 0;
subject to v5_6_7:
	-z5_6_6 - y5_5_7 - x4_6_7 + x5_6_7 + y5_6_7 + z5_6_7 - 1.0 = 0;
subject to v6_6_7:
	-z6_6_6 - y6_5_7 - x5_6_7 + x6_6_7 + y6_6_7 + z6_6_7 - 1.0 = 0;
subject to v7_6_7:
	-z7_6_6 - y7_5_7 - x6_6_7 + x7_6_7 + y7_6_7 + z7_6_7 - 1.0 = 0;
subject to v8_6_7:
	-z8_6_6 - y8_5_7 - x7_6_7 + x8_6_7 + y8_6_7 + z8_6_7 - 1.0 = 0;
subject to v9_6_7:
	-z9_6_6 - y9_5_7 - x8_6_7 + x9_6_7 + y9_6_7 + z9_6_7 - 1.0 = 0;
subject to v10_6_7:
	-x9_6_7 - z10_6_6 - y10_5_7 + y10_6_7 + z10_6_7 + y11_6_7 + z11_6_7 - 1.0 = 0;
subject to v1_7_7:
	-z1_7_6 - y1_6_7 + x1_7_7 + y1_7_7 + z1_7_7 + y0_7_7 + z0_7_7 - 1.0 = 0;
subject to v2_7_7:
	-z2_7_6 - y2_6_7 - x1_7_7 + x2_7_7 + y2_7_7 + z2_7_7 - 1.0 = 0;
subject to v3_7_7:
	-z3_7_6 - y3_6_7 - x2_7_7 + x3_7_7 + y3_7_7 + z3_7_7 - 1.0 = 0;
subject to v4_7_7:
	-z4_7_6 - y4_6_7 - x3_7_7 + x4_7_7 + y4_7_7 + z4_7_7 - 1.0 = 0;
subject to v5_7_7:
	-z5_7_6 - y5_6_7 - x4_7_7 + x5_7_7 + y5_7_7 + z5_7_7 - 1.0 = 0;
subject to v6_7_7:
	-z6_7_6 - y6_6_7 - x5_7_7 + x6_7_7 + y6_7_7 + z6_7_7 - 1.0 = 0;
subject to v7_7_7:
	-z7_7_6 - y7_6_7 - x6_7_7 + x7_7_7 + y7_7_7 + z7_7_7 - 1.0 = 0;
subject to v8_7_7:
	-z8_7_6 - y8_6_7 - x7_7_7 + x8_7_7 + y8_7_7 + z8_7_7 - 1.0 = 0;
subject to v9_7_7:
	-z9_7_6 - y9_6_7 - x8_7_7 + x9_7_7 + y9_7_7 + z9_7_7 - 1.0 = 0;
subject to v10_7_7:
	-x9_7_7 - z10_7_6 - y10_6_7 + y10_7_7 + z10_7_7 + y11_7_7 + z11_7_7 - 1.0 = 0;
subject to v1_8_7:
	-z1_8_6 - y1_7_7 + x1_8_7 + y1_8_7 + z1_8_7 + y0_8_7 + z0_8_7 - 1.0 = 0;
subject to v2_8_7:
	-z2_8_6 - y2_7_7 - x1_8_7 + x2_8_7 + y2_8_7 + z2_8_7 - 1.0 = 0;
subject to v3_8_7:
	-z3_8_6 - y3_7_7 - x2_8_7 + x3_8_7 + y3_8_7 + z3_8_7 - 1.0 = 0;
subject to v4_8_7:
	-z4_8_6 - y4_7_7 - x3_8_7 + x4_8_7 + y4_8_7 + z4_8_7 - 1.0 = 0;
subject to v5_8_7:
	-z5_8_6 - y5_7_7 - x4_8_7 + x5_8_7 + y5_8_7 + z5_8_7 - 1.0 = 0;
subject to v6_8_7:
	-z6_8_6 - y6_7_7 - x5_8_7 + x6_8_7 + y6_8_7 + z6_8_7 - 1.0 = 0;
subject to v7_8_7:
	-z7_8_6 - y7_7_7 - x6_8_7 + x7_8_7 + y7_8_7 + z7_8_7 - 1.0 = 0;
subject to v8_8_7:
	-z8_8_6 - y8_7_7 - x7_8_7 + x8_8_7 + y8_8_7 + z8_8_7 - 1.0 = 0;
subject to v9_8_7:
	-z9_8_6 - y9_7_7 - x8_8_7 + x9_8_7 + y9_8_7 + z9_8_7 - 1.0 = 0;
subject to v10_8_7:
	-x9_8_7 - z10_8_6 - y10_7_7 + y10_8_7 + z10_8_7 + y11_8_7 + z11_8_7 - 1.0 = 0;
subject to v1_9_7:
	-z1_9_6 - y1_8_7 + x1_9_7 + y1_9_7 + z1_9_7 + y0_9_7 + z0_9_7 - 1.0 = 0;
subject to v2_9_7:
	-z2_9_6 - y2_8_7 - x1_9_7 + x2_9_7 + y2_9_7 + z2_9_7 - 1.0 = 0;
subject to v3_9_7:
	-z3_9_6 - y3_8_7 - x2_9_7 + x3_9_7 + y3_9_7 + z3_9_7 - 1.0 = 0;
subject to v4_9_7:
	-z4_9_6 - y4_8_7 - x3_9_7 + x4_9_7 + y4_9_7 + z4_9_7 - 1.0 = 0;
subject to v5_9_7:
	-z5_9_6 - y5_8_7 - x4_9_7 + x5_9_7 + y5_9_7 + z5_9_7 - 1.0 = 0;
subject to v6_9_7:
	-z6_9_6 - y6_8_7 - x5_9_7 + x6_9_7 + y6_9_7 + z6_9_7 - 1.0 = 0;
subject to v7_9_7:
	-z7_9_6 - y7_8_7 - x6_9_7 + x7_9_7 + y7_9_7 + z7_9_7 - 1.0 = 0;
subject to v8_9_7:
	-z8_9_6 - y8_8_7 - x7_9_7 + x8_9_7 + y8_9_7 + z8_9_7 - 1.0 = 0;
subject to v9_9_7:
	-z9_9_6 - y9_8_7 - x8_9_7 + x9_9_7 + y9_9_7 + z9_9_7 - 1.0 = 0;
subject to v10_9_7:
	-x9_9_7 - z10_9_6 - y10_8_7 + y10_9_7 + z10_9_7 + y11_9_7 + z11_9_7 - 1.0 = 0;
subject to v1_10_7:
	-y1_9_7 - z1_10_6 + x1_10_7 + z1_10_7 + y0_10_7 + z0_10_7 + x1_11_7 + z1_11_7 - 
	1.0 = 0;
subject to v2_10_7:
	-y2_9_7 - z2_10_6 - x1_10_7 + x2_10_7 + z2_10_7 + x2_11_7 + z2_11_7 - 1.0 = 0;
subject to v3_10_7:
	-y3_9_7 - z3_10_6 - x2_10_7 + x3_10_7 + z3_10_7 + x3_11_7 + z3_11_7 - 1.0 = 0;
subject to v4_10_7:
	-y4_9_7 - z4_10_6 - x3_10_7 + x4_10_7 + z4_10_7 + x4_11_7 + z4_11_7 - 1.0 = 0;
subject to v5_10_7:
	-y5_9_7 - z5_10_6 - x4_10_7 + x5_10_7 + z5_10_7 + x5_11_7 + z5_11_7 - 1.0 = 0;
subject to v6_10_7:
	-y6_9_7 - z6_10_6 - x5_10_7 + x6_10_7 + z6_10_7 + x6_11_7 + z6_11_7 - 1.0 = 0;
subject to v7_10_7:
	-y7_9_7 - z7_10_6 - x6_10_7 + x7_10_7 + z7_10_7 + x7_11_7 + z7_11_7 - 1.0 = 0;
subject to v8_10_7:
	-y8_9_7 - z8_10_6 - x7_10_7 + x8_10_7 + z8_10_7 + x8_11_7 + z8_11_7 - 1.0 = 0;
subject to v9_10_7:
	-y9_9_7 - z9_10_6 - x8_10_7 + x9_10_7 + z9_10_7 + x9_11_7 + z9_11_7 - 1.0 = 0;
subject to v10_10_7:
	-y10_9_7 - x9_10_7 + y11_10_7 + z11_10_7 + x10_11_7 + z10_11_7 - 1.0 = 0;
subject to v1_1_8:
	-z1_1_7 + x1_1_8 + y1_1_8 + z1_1_8 + y0_1_8 + z0_1_8 + x1_0_8 + z1_0_8 - 1.0 = 
	0;
subject to v2_1_8:
	-z2_1_7 - x1_1_8 + x2_1_8 + y2_1_8 + z2_1_8 + x2_0_8 + z2_0_8 - 1.0 = 0;
subject to v3_1_8:
	-z3_1_7 - x2_1_8 + x3_1_8 + y3_1_8 + z3_1_8 + x3_0_8 + z3_0_8 - 1.0 = 0;
subject to v4_1_8:
	-z4_1_7 - x3_1_8 + x4_1_8 + y4_1_8 + z4_1_8 + x4_0_8 + z4_0_8 - 1.0 = 0;
subject to v5_1_8:
	-z5_1_7 - x4_1_8 + x5_1_8 + y5_1_8 + z5_1_8 + x5_0_8 + z5_0_8 - 1.0 = 0;
subject to v6_1_8:
	-z6_1_7 - x5_1_8 + x6_1_8 + y6_1_8 + z6_1_8 + x6_0_8 + z6_0_8 - 1.0 = 0;
subject to v7_1_8:
	-z7_1_7 - x6_1_8 + x7_1_8 + y7_1_8 + z7_1_8 + x7_0_8 + z7_0_8 - 1.0 = 0;
subject to v8_1_8:
	-z8_1_7 - x7_1_8 + x8_1_8 + y8_1_8 + z8_1_8 + x8_0_8 + z8_0_8 - 1.0 = 0;
subject to v9_1_8:
	-z9_1_7 - x8_1_8 + x9_1_8 + y9_1_8 + z9_1_8 + x9_0_8 + z9_0_8 - 1.0 = 0;
subject to v10_1_8:
	-x9_1_8 - z10_1_7 + y10_1_8 + z10_1_8 + y11_1_8 + z11_1_8 + x10_0_8 + z10_0_8 - 
	1.0 = 0;
subject to v1_2_8:
	-z1_2_7 - y1_1_8 + x1_2_8 + y1_2_8 + z1_2_8 + y0_2_8 + z0_2_8 - 1.0 = 0;
subject to v2_2_8:
	-z2_2_7 - y2_1_8 - x1_2_8 + x2_2_8 + y2_2_8 + z2_2_8 - 1.0 = 0;
subject to v3_2_8:
	-z3_2_7 - y3_1_8 - x2_2_8 + x3_2_8 + y3_2_8 + z3_2_8 - 1.0 = 0;
subject to v4_2_8:
	-z4_2_7 - y4_1_8 - x3_2_8 + x4_2_8 + y4_2_8 + z4_2_8 - 1.0 = 0;
subject to v5_2_8:
	-z5_2_7 - y5_1_8 - x4_2_8 + x5_2_8 + y5_2_8 + z5_2_8 - 1.0 = 0;
subject to v6_2_8:
	-z6_2_7 - y6_1_8 - x5_2_8 + x6_2_8 + y6_2_8 + z6_2_8 - 1.0 = 0;
subject to v7_2_8:
	-z7_2_7 - y7_1_8 - x6_2_8 + x7_2_8 + y7_2_8 + z7_2_8 - 1.0 = 0;
subject to v8_2_8:
	-z8_2_7 - y8_1_8 - x7_2_8 + x8_2_8 + y8_2_8 + z8_2_8 - 1.0 = 0;
subject to v9_2_8:
	-z9_2_7 - y9_1_8 - x8_2_8 + x9_2_8 + y9_2_8 + z9_2_8 - 1.0 = 0;
subject to v10_2_8:
	-x9_2_8 - z10_2_7 - y10_1_8 + y10_2_8 + z10_2_8 + y11_2_8 + z11_2_8 - 1.0 = 0;
subject to v1_3_8:
	-z1_3_7 - y1_2_8 + x1_3_8 + y1_3_8 + z1_3_8 + y0_3_8 + z0_3_8 - 1.0 = 0;
subject to v2_3_8:
	-z2_3_7 - y2_2_8 - x1_3_8 + x2_3_8 + y2_3_8 + z2_3_8 - 1.0 = 0;
subject to v3_3_8:
	-z3_3_7 - y3_2_8 - x2_3_8 + x3_3_8 + y3_3_8 + z3_3_8 - 1.0 = 0;
subject to v4_3_8:
	-z4_3_7 - y4_2_8 - x3_3_8 + x4_3_8 + y4_3_8 + z4_3_8 - 1.0 = 0;
subject to v5_3_8:
	-z5_3_7 - y5_2_8 - x4_3_8 + x5_3_8 + y5_3_8 + z5_3_8 - 1.0 = 0;
subject to v6_3_8:
	-z6_3_7 - y6_2_8 - x5_3_8 + x6_3_8 + y6_3_8 + z6_3_8 - 1.0 = 0;
subject to v7_3_8:
	-z7_3_7 - y7_2_8 - x6_3_8 + x7_3_8 + y7_3_8 + z7_3_8 - 1.0 = 0;
subject to v8_3_8:
	-z8_3_7 - y8_2_8 - x7_3_8 + x8_3_8 + y8_3_8 + z8_3_8 - 1.0 = 0;
subject to v9_3_8:
	-z9_3_7 - y9_2_8 - x8_3_8 + x9_3_8 + y9_3_8 + z9_3_8 - 1.0 = 0;
subject to v10_3_8:
	-x9_3_8 - z10_3_7 - y10_2_8 + y10_3_8 + z10_3_8 + y11_3_8 + z11_3_8 - 1.0 = 0;
subject to v1_4_8:
	-z1_4_7 - y1_3_8 + x1_4_8 + y1_4_8 + z1_4_8 + y0_4_8 + z0_4_8 - 1.0 = 0;
subject to v2_4_8:
	-z2_4_7 - y2_3_8 - x1_4_8 + x2_4_8 + y2_4_8 + z2_4_8 - 1.0 = 0;
subject to v3_4_8:
	-z3_4_7 - y3_3_8 - x2_4_8 + x3_4_8 + y3_4_8 + z3_4_8 - 1.0 = 0;
subject to v4_4_8:
	-z4_4_7 - y4_3_8 - x3_4_8 + x4_4_8 + y4_4_8 + z4_4_8 - 1.0 = 0;
subject to v5_4_8:
	-z5_4_7 - y5_3_8 - x4_4_8 + x5_4_8 + y5_4_8 + z5_4_8 - 1.0 = 0;
subject to v6_4_8:
	-z6_4_7 - y6_3_8 - x5_4_8 + x6_4_8 + y6_4_8 + z6_4_8 - 1.0 = 0;
subject to v7_4_8:
	-z7_4_7 - y7_3_8 - x6_4_8 + x7_4_8 + y7_4_8 + z7_4_8 - 1.0 = 0;
subject to v8_4_8:
	-z8_4_7 - y8_3_8 - x7_4_8 + x8_4_8 + y8_4_8 + z8_4_8 - 1.0 = 0;
subject to v9_4_8:
	-z9_4_7 - y9_3_8 - x8_4_8 + x9_4_8 + y9_4_8 + z9_4_8 - 1.0 = 0;
subject to v10_4_8:
	-x9_4_8 - z10_4_7 - y10_3_8 + y10_4_8 + z10_4_8 + y11_4_8 + z11_4_8 - 1.0 = 0;
subject to v1_5_8:
	-z1_5_7 - y1_4_8 + x1_5_8 + y1_5_8 + z1_5_8 + y0_5_8 + z0_5_8 - 1.0 = 0;
subject to v2_5_8:
	-z2_5_7 - y2_4_8 - x1_5_8 + x2_5_8 + y2_5_8 + z2_5_8 - 1.0 = 0;
subject to v3_5_8:
	-z3_5_7 - y3_4_8 - x2_5_8 + x3_5_8 + y3_5_8 + z3_5_8 - 1.0 = 0;
subject to v4_5_8:
	-z4_5_7 - y4_4_8 - x3_5_8 + x4_5_8 + y4_5_8 + z4_5_8 - 1.0 = 0;
subject to v5_5_8:
	-z5_5_7 - y5_4_8 - x4_5_8 + x5_5_8 + y5_5_8 + z5_5_8 - 1.0 = 0;
subject to v6_5_8:
	-z6_5_7 - y6_4_8 - x5_5_8 + x6_5_8 + y6_5_8 + z6_5_8 - 1.0 = 0;
subject to v7_5_8:
	-z7_5_7 - y7_4_8 - x6_5_8 + x7_5_8 + y7_5_8 + z7_5_8 - 1.0 = 0;
subject to v8_5_8:
	-z8_5_7 - y8_4_8 - x7_5_8 + x8_5_8 + y8_5_8 + z8_5_8 - 1.0 = 0;
subject to v9_5_8:
	-z9_5_7 - y9_4_8 - x8_5_8 + x9_5_8 + y9_5_8 + z9_5_8 - 1.0 = 0;
subject to v10_5_8:
	-x9_5_8 - z10_5_7 - y10_4_8 + y10_5_8 + z10_5_8 + y11_5_8 + z11_5_8 - 1.0 = 0;
subject to v1_6_8:
	-z1_6_7 - y1_5_8 + x1_6_8 + y1_6_8 + z1_6_8 + y0_6_8 + z0_6_8 - 1.0 = 0;
subject to v2_6_8:
	-z2_6_7 - y2_5_8 - x1_6_8 + x2_6_8 + y2_6_8 + z2_6_8 - 1.0 = 0;
subject to v3_6_8:
	-z3_6_7 - y3_5_8 - x2_6_8 + x3_6_8 + y3_6_8 + z3_6_8 - 1.0 = 0;
subject to v4_6_8:
	-z4_6_7 - y4_5_8 - x3_6_8 + x4_6_8 + y4_6_8 + z4_6_8 - 1.0 = 0;
subject to v5_6_8:
	-z5_6_7 - y5_5_8 - x4_6_8 + x5_6_8 + y5_6_8 + z5_6_8 - 1.0 = 0;
subject to v6_6_8:
	-z6_6_7 - y6_5_8 - x5_6_8 + x6_6_8 + y6_6_8 + z6_6_8 - 1.0 = 0;
subject to v7_6_8:
	-z7_6_7 - y7_5_8 - x6_6_8 + x7_6_8 + y7_6_8 + z7_6_8 - 1.0 = 0;
subject to v8_6_8:
	-z8_6_7 - y8_5_8 - x7_6_8 + x8_6_8 + y8_6_8 + z8_6_8 - 1.0 = 0;
subject to v9_6_8:
	-z9_6_7 - y9_5_8 - x8_6_8 + x9_6_8 + y9_6_8 + z9_6_8 - 1.0 = 0;
subject to v10_6_8:
	-x9_6_8 - z10_6_7 - y10_5_8 + y10_6_8 + z10_6_8 + y11_6_8 + z11_6_8 - 1.0 = 0;
subject to v1_7_8:
	-z1_7_7 - y1_6_8 + x1_7_8 + y1_7_8 + z1_7_8 + y0_7_8 + z0_7_8 - 1.0 = 0;
subject to v2_7_8:
	-z2_7_7 - y2_6_8 - x1_7_8 + x2_7_8 + y2_7_8 + z2_7_8 - 1.0 = 0;
subject to v3_7_8:
	-z3_7_7 - y3_6_8 - x2_7_8 + x3_7_8 + y3_7_8 + z3_7_8 - 1.0 = 0;
subject to v4_7_8:
	-z4_7_7 - y4_6_8 - x3_7_8 + x4_7_8 + y4_7_8 + z4_7_8 - 1.0 = 0;
subject to v5_7_8:
	-z5_7_7 - y5_6_8 - x4_7_8 + x5_7_8 + y5_7_8 + z5_7_8 - 1.0 = 0;
subject to v6_7_8:
	-z6_7_7 - y6_6_8 - x5_7_8 + x6_7_8 + y6_7_8 + z6_7_8 - 1.0 = 0;
subject to v7_7_8:
	-z7_7_7 - y7_6_8 - x6_7_8 + x7_7_8 + y7_7_8 + z7_7_8 - 1.0 = 0;
subject to v8_7_8:
	-z8_7_7 - y8_6_8 - x7_7_8 + x8_7_8 + y8_7_8 + z8_7_8 - 1.0 = 0;
subject to v9_7_8:
	-z9_7_7 - y9_6_8 - x8_7_8 + x9_7_8 + y9_7_8 + z9_7_8 - 1.0 = 0;
subject to v10_7_8:
	-x9_7_8 - z10_7_7 - y10_6_8 + y10_7_8 + z10_7_8 + y11_7_8 + z11_7_8 - 1.0 = 0;
subject to v1_8_8:
	-z1_8_7 - y1_7_8 + x1_8_8 + y1_8_8 + z1_8_8 + y0_8_8 + z0_8_8 - 1.0 = 0;
subject to v2_8_8:
	-z2_8_7 - y2_7_8 - x1_8_8 + x2_8_8 + y2_8_8 + z2_8_8 - 1.0 = 0;
subject to v3_8_8:
	-z3_8_7 - y3_7_8 - x2_8_8 + x3_8_8 + y3_8_8 + z3_8_8 - 1.0 = 0;
subject to v4_8_8:
	-z4_8_7 - y4_7_8 - x3_8_8 + x4_8_8 + y4_8_8 + z4_8_8 - 1.0 = 0;
subject to v5_8_8:
	-z5_8_7 - y5_7_8 - x4_8_8 + x5_8_8 + y5_8_8 + z5_8_8 - 1.0 = 0;
subject to v6_8_8:
	-z6_8_7 - y6_7_8 - x5_8_8 + x6_8_8 + y6_8_8 + z6_8_8 - 1.0 = 0;
subject to v7_8_8:
	-z7_8_7 - y7_7_8 - x6_8_8 + x7_8_8 + y7_8_8 + z7_8_8 - 1.0 = 0;
subject to v8_8_8:
	-z8_8_7 - y8_7_8 - x7_8_8 + x8_8_8 + y8_8_8 + z8_8_8 - 1.0 = 0;
subject to v9_8_8:
	-z9_8_7 - y9_7_8 - x8_8_8 + x9_8_8 + y9_8_8 + z9_8_8 - 1.0 = 0;
subject to v10_8_8:
	-x9_8_8 - z10_8_7 - y10_7_8 + y10_8_8 + z10_8_8 + y11_8_8 + z11_8_8 - 1.0 = 0;
subject to v1_9_8:
	-z1_9_7 - y1_8_8 + x1_9_8 + y1_9_8 + z1_9_8 + y0_9_8 + z0_9_8 - 1.0 = 0;
subject to v2_9_8:
	-z2_9_7 - y2_8_8 - x1_9_8 + x2_9_8 + y2_9_8 + z2_9_8 - 1.0 = 0;
subject to v3_9_8:
	-z3_9_7 - y3_8_8 - x2_9_8 + x3_9_8 + y3_9_8 + z3_9_8 - 1.0 = 0;
subject to v4_9_8:
	-z4_9_7 - y4_8_8 - x3_9_8 + x4_9_8 + y4_9_8 + z4_9_8 - 1.0 = 0;
subject to v5_9_8:
	-z5_9_7 - y5_8_8 - x4_9_8 + x5_9_8 + y5_9_8 + z5_9_8 - 1.0 = 0;
subject to v6_9_8:
	-z6_9_7 - y6_8_8 - x5_9_8 + x6_9_8 + y6_9_8 + z6_9_8 - 1.0 = 0;
subject to v7_9_8:
	-z7_9_7 - y7_8_8 - x6_9_8 + x7_9_8 + y7_9_8 + z7_9_8 - 1.0 = 0;
subject to v8_9_8:
	-z8_9_7 - y8_8_8 - x7_9_8 + x8_9_8 + y8_9_8 + z8_9_8 - 1.0 = 0;
subject to v9_9_8:
	-z9_9_7 - y9_8_8 - x8_9_8 + x9_9_8 + y9_9_8 + z9_9_8 - 1.0 = 0;
subject to v10_9_8:
	-x9_9_8 - z10_9_7 - y10_8_8 + y10_9_8 + z10_9_8 + y11_9_8 + z11_9_8 - 1.0 = 0;
subject to v1_10_8:
	-y1_9_8 - z1_10_7 + x1_10_8 + z1_10_8 + y0_10_8 + z0_10_8 + x1_11_8 + z1_11_8 - 
	1.0 = 0;
subject to v2_10_8:
	-y2_9_8 - z2_10_7 - x1_10_8 + x2_10_8 + z2_10_8 + x2_11_8 + z2_11_8 - 1.0 = 0;
subject to v3_10_8:
	-y3_9_8 - z3_10_7 - x2_10_8 + x3_10_8 + z3_10_8 + x3_11_8 + z3_11_8 - 1.0 = 0;
subject to v4_10_8:
	-y4_9_8 - z4_10_7 - x3_10_8 + x4_10_8 + z4_10_8 + x4_11_8 + z4_11_8 - 1.0 = 0;
subject to v5_10_8:
	-y5_9_8 - z5_10_7 - x4_10_8 + x5_10_8 + z5_10_8 + x5_11_8 + z5_11_8 - 1.0 = 0;
subject to v6_10_8:
	-y6_9_8 - z6_10_7 - x5_10_8 + x6_10_8 + z6_10_8 + x6_11_8 + z6_11_8 - 1.0 = 0;
subject to v7_10_8:
	-y7_9_8 - z7_10_7 - x6_10_8 + x7_10_8 + z7_10_8 + x7_11_8 + z7_11_8 - 1.0 = 0;
subject to v8_10_8:
	-y8_9_8 - z8_10_7 - x7_10_8 + x8_10_8 + z8_10_8 + x8_11_8 + z8_11_8 - 1.0 = 0;
subject to v9_10_8:
	-y9_9_8 - z9_10_7 - x8_10_8 + x9_10_8 + z9_10_8 + x9_11_8 + z9_11_8 - 1.0 = 0;
subject to v10_10_8:
	-y10_9_8 - x9_10_8 + y11_10_8 + z11_10_8 + x10_11_8 + z10_11_8 - 1.0 = 0;
subject to v1_1_9:
	-z1_1_8 + x1_1_9 + y1_1_9 + z1_1_9 + y0_1_9 + z0_1_9 + x1_0_9 + z1_0_9 - 1.0 = 
	0;
subject to v2_1_9:
	-z2_1_8 - x1_1_9 + x2_1_9 + y2_1_9 + z2_1_9 + x2_0_9 + z2_0_9 - 1.0 = 0;
subject to v3_1_9:
	-z3_1_8 - x2_1_9 + x3_1_9 + y3_1_9 + z3_1_9 + x3_0_9 + z3_0_9 - 1.0 = 0;
subject to v4_1_9:
	-z4_1_8 - x3_1_9 + x4_1_9 + y4_1_9 + z4_1_9 + x4_0_9 + z4_0_9 - 1.0 = 0;
subject to v5_1_9:
	-z5_1_8 - x4_1_9 + x5_1_9 + y5_1_9 + z5_1_9 + x5_0_9 + z5_0_9 - 1.0 = 0;
subject to v6_1_9:
	-z6_1_8 - x5_1_9 + x6_1_9 + y6_1_9 + z6_1_9 + x6_0_9 + z6_0_9 - 1.0 = 0;
subject to v7_1_9:
	-z7_1_8 - x6_1_9 + x7_1_9 + y7_1_9 + z7_1_9 + x7_0_9 + z7_0_9 - 1.0 = 0;
subject to v8_1_9:
	-z8_1_8 - x7_1_9 + x8_1_9 + y8_1_9 + z8_1_9 + x8_0_9 + z8_0_9 - 1.0 = 0;
subject to v9_1_9:
	-z9_1_8 - x8_1_9 + x9_1_9 + y9_1_9 + z9_1_9 + x9_0_9 + z9_0_9 - 1.0 = 0;
subject to v10_1_9:
	-x9_1_9 - z10_1_8 + y10_1_9 + z10_1_9 + y11_1_9 + z11_1_9 + x10_0_9 + z10_0_9 - 
	1.0 = 0;
subject to v1_2_9:
	-z1_2_8 - y1_1_9 + x1_2_9 + y1_2_9 + z1_2_9 + y0_2_9 + z0_2_9 - 1.0 = 0;
subject to v2_2_9:
	-z2_2_8 - y2_1_9 - x1_2_9 + x2_2_9 + y2_2_9 + z2_2_9 - 1.0 = 0;
subject to v3_2_9:
	-z3_2_8 - y3_1_9 - x2_2_9 + x3_2_9 + y3_2_9 + z3_2_9 - 1.0 = 0;
subject to v4_2_9:
	-z4_2_8 - y4_1_9 - x3_2_9 + x4_2_9 + y4_2_9 + z4_2_9 - 1.0 = 0;
subject to v5_2_9:
	-z5_2_8 - y5_1_9 - x4_2_9 + x5_2_9 + y5_2_9 + z5_2_9 - 1.0 = 0;
subject to v6_2_9:
	-z6_2_8 - y6_1_9 - x5_2_9 + x6_2_9 + y6_2_9 + z6_2_9 - 1.0 = 0;
subject to v7_2_9:
	-z7_2_8 - y7_1_9 - x6_2_9 + x7_2_9 + y7_2_9 + z7_2_9 - 1.0 = 0;
subject to v8_2_9:
	-z8_2_8 - y8_1_9 - x7_2_9 + x8_2_9 + y8_2_9 + z8_2_9 - 1.0 = 0;
subject to v9_2_9:
	-z9_2_8 - y9_1_9 - x8_2_9 + x9_2_9 + y9_2_9 + z9_2_9 - 1.0 = 0;
subject to v10_2_9:
	-x9_2_9 - z10_2_8 - y10_1_9 + y10_2_9 + z10_2_9 + y11_2_9 + z11_2_9 - 1.0 = 0;
subject to v1_3_9:
	-z1_3_8 - y1_2_9 + x1_3_9 + y1_3_9 + z1_3_9 + y0_3_9 + z0_3_9 - 1.0 = 0;
subject to v2_3_9:
	-z2_3_8 - y2_2_9 - x1_3_9 + x2_3_9 + y2_3_9 + z2_3_9 - 1.0 = 0;
subject to v3_3_9:
	-z3_3_8 - y3_2_9 - x2_3_9 + x3_3_9 + y3_3_9 + z3_3_9 - 1.0 = 0;
subject to v4_3_9:
	-z4_3_8 - y4_2_9 - x3_3_9 + x4_3_9 + y4_3_9 + z4_3_9 - 1.0 = 0;
subject to v5_3_9:
	-z5_3_8 - y5_2_9 - x4_3_9 + x5_3_9 + y5_3_9 + z5_3_9 - 1.0 = 0;
subject to v6_3_9:
	-z6_3_8 - y6_2_9 - x5_3_9 + x6_3_9 + y6_3_9 + z6_3_9 - 1.0 = 0;
subject to v7_3_9:
	-z7_3_8 - y7_2_9 - x6_3_9 + x7_3_9 + y7_3_9 + z7_3_9 - 1.0 = 0;
subject to v8_3_9:
	-z8_3_8 - y8_2_9 - x7_3_9 + x8_3_9 + y8_3_9 + z8_3_9 - 1.0 = 0;
subject to v9_3_9:
	-z9_3_8 - y9_2_9 - x8_3_9 + x9_3_9 + y9_3_9 + z9_3_9 - 1.0 = 0;
subject to v10_3_9:
	-x9_3_9 - z10_3_8 - y10_2_9 + y10_3_9 + z10_3_9 + y11_3_9 + z11_3_9 - 1.0 = 0;
subject to v1_4_9:
	-z1_4_8 - y1_3_9 + x1_4_9 + y1_4_9 + z1_4_9 + y0_4_9 + z0_4_9 - 1.0 = 0;
subject to v2_4_9:
	-z2_4_8 - y2_3_9 - x1_4_9 + x2_4_9 + y2_4_9 + z2_4_9 - 1.0 = 0;
subject to v3_4_9:
	-z3_4_8 - y3_3_9 - x2_4_9 + x3_4_9 + y3_4_9 + z3_4_9 - 1.0 = 0;
subject to v4_4_9:
	-z4_4_8 - y4_3_9 - x3_4_9 + x4_4_9 + y4_4_9 + z4_4_9 - 1.0 = 0;
subject to v5_4_9:
	-z5_4_8 - y5_3_9 - x4_4_9 + x5_4_9 + y5_4_9 + z5_4_9 - 1.0 = 0;
subject to v6_4_9:
	-z6_4_8 - y6_3_9 - x5_4_9 + x6_4_9 + y6_4_9 + z6_4_9 - 1.0 = 0;
subject to v7_4_9:
	-z7_4_8 - y7_3_9 - x6_4_9 + x7_4_9 + y7_4_9 + z7_4_9 - 1.0 = 0;
subject to v8_4_9:
	-z8_4_8 - y8_3_9 - x7_4_9 + x8_4_9 + y8_4_9 + z8_4_9 - 1.0 = 0;
subject to v9_4_9:
	-z9_4_8 - y9_3_9 - x8_4_9 + x9_4_9 + y9_4_9 + z9_4_9 - 1.0 = 0;
subject to v10_4_9:
	-x9_4_9 - z10_4_8 - y10_3_9 + y10_4_9 + z10_4_9 + y11_4_9 + z11_4_9 - 1.0 = 0;
subject to v1_5_9:
	-z1_5_8 - y1_4_9 + x1_5_9 + y1_5_9 + z1_5_9 + y0_5_9 + z0_5_9 - 1.0 = 0;
subject to v2_5_9:
	-z2_5_8 - y2_4_9 - x1_5_9 + x2_5_9 + y2_5_9 + z2_5_9 - 1.0 = 0;
subject to v3_5_9:
	-z3_5_8 - y3_4_9 - x2_5_9 + x3_5_9 + y3_5_9 + z3_5_9 - 1.0 = 0;
subject to v4_5_9:
	-z4_5_8 - y4_4_9 - x3_5_9 + x4_5_9 + y4_5_9 + z4_5_9 - 1.0 = 0;
subject to v5_5_9:
	-z5_5_8 - y5_4_9 - x4_5_9 + x5_5_9 + y5_5_9 + z5_5_9 - 1.0 = 0;
subject to v6_5_9:
	-z6_5_8 - y6_4_9 - x5_5_9 + x6_5_9 + y6_5_9 + z6_5_9 - 1.0 = 0;
subject to v7_5_9:
	-z7_5_8 - y7_4_9 - x6_5_9 + x7_5_9 + y7_5_9 + z7_5_9 - 1.0 = 0;
subject to v8_5_9:
	-z8_5_8 - y8_4_9 - x7_5_9 + x8_5_9 + y8_5_9 + z8_5_9 - 1.0 = 0;
subject to v9_5_9:
	-z9_5_8 - y9_4_9 - x8_5_9 + x9_5_9 + y9_5_9 + z9_5_9 - 1.0 = 0;
subject to v10_5_9:
	-x9_5_9 - z10_5_8 - y10_4_9 + y10_5_9 + z10_5_9 + y11_5_9 + z11_5_9 - 1.0 = 0;
subject to v1_6_9:
	-z1_6_8 - y1_5_9 + x1_6_9 + y1_6_9 + z1_6_9 + y0_6_9 + z0_6_9 - 1.0 = 0;
subject to v2_6_9:
	-z2_6_8 - y2_5_9 - x1_6_9 + x2_6_9 + y2_6_9 + z2_6_9 - 1.0 = 0;
subject to v3_6_9:
	-z3_6_8 - y3_5_9 - x2_6_9 + x3_6_9 + y3_6_9 + z3_6_9 - 1.0 = 0;
subject to v4_6_9:
	-z4_6_8 - y4_5_9 - x3_6_9 + x4_6_9 + y4_6_9 + z4_6_9 - 1.0 = 0;
subject to v5_6_9:
	-z5_6_8 - y5_5_9 - x4_6_9 + x5_6_9 + y5_6_9 + z5_6_9 - 1.0 = 0;
subject to v6_6_9:
	-z6_6_8 - y6_5_9 - x5_6_9 + x6_6_9 + y6_6_9 + z6_6_9 - 1.0 = 0;
subject to v7_6_9:
	-z7_6_8 - y7_5_9 - x6_6_9 + x7_6_9 + y7_6_9 + z7_6_9 - 1.0 = 0;
subject to v8_6_9:
	-z8_6_8 - y8_5_9 - x7_6_9 + x8_6_9 + y8_6_9 + z8_6_9 - 1.0 = 0;
subject to v9_6_9:
	-z9_6_8 - y9_5_9 - x8_6_9 + x9_6_9 + y9_6_9 + z9_6_9 - 1.0 = 0;
subject to v10_6_9:
	-x9_6_9 - z10_6_8 - y10_5_9 + y10_6_9 + z10_6_9 + y11_6_9 + z11_6_9 - 1.0 = 0;
subject to v1_7_9:
	-z1_7_8 - y1_6_9 + x1_7_9 + y1_7_9 + z1_7_9 + y0_7_9 + z0_7_9 - 1.0 = 0;
subject to v2_7_9:
	-z2_7_8 - y2_6_9 - x1_7_9 + x2_7_9 + y2_7_9 + z2_7_9 - 1.0 = 0;
subject to v3_7_9:
	-z3_7_8 - y3_6_9 - x2_7_9 + x3_7_9 + y3_7_9 + z3_7_9 - 1.0 = 0;
subject to v4_7_9:
	-z4_7_8 - y4_6_9 - x3_7_9 + x4_7_9 + y4_7_9 + z4_7_9 - 1.0 = 0;
subject to v5_7_9:
	-z5_7_8 - y5_6_9 - x4_7_9 + x5_7_9 + y5_7_9 + z5_7_9 - 1.0 = 0;
subject to v6_7_9:
	-z6_7_8 - y6_6_9 - x5_7_9 + x6_7_9 + y6_7_9 + z6_7_9 - 1.0 = 0;
subject to v7_7_9:
	-z7_7_8 - y7_6_9 - x6_7_9 + x7_7_9 + y7_7_9 + z7_7_9 - 1.0 = 0;
subject to v8_7_9:
	-z8_7_8 - y8_6_9 - x7_7_9 + x8_7_9 + y8_7_9 + z8_7_9 - 1.0 = 0;
subject to v9_7_9:
	-z9_7_8 - y9_6_9 - x8_7_9 + x9_7_9 + y9_7_9 + z9_7_9 - 1.0 = 0;
subject to v10_7_9:
	-x9_7_9 - z10_7_8 - y10_6_9 + y10_7_9 + z10_7_9 + y11_7_9 + z11_7_9 - 1.0 = 0;
subject to v1_8_9:
	-z1_8_8 - y1_7_9 + x1_8_9 + y1_8_9 + z1_8_9 + y0_8_9 + z0_8_9 - 1.0 = 0;
subject to v2_8_9:
	-z2_8_8 - y2_7_9 - x1_8_9 + x2_8_9 + y2_8_9 + z2_8_9 - 1.0 = 0;
subject to v3_8_9:
	-z3_8_8 - y3_7_9 - x2_8_9 + x3_8_9 + y3_8_9 + z3_8_9 - 1.0 = 0;
subject to v4_8_9:
	-z4_8_8 - y4_7_9 - x3_8_9 + x4_8_9 + y4_8_9 + z4_8_9 - 1.0 = 0;
subject to v5_8_9:
	-z5_8_8 - y5_7_9 - x4_8_9 + x5_8_9 + y5_8_9 + z5_8_9 - 1.0 = 0;
subject to v6_8_9:
	-z6_8_8 - y6_7_9 - x5_8_9 + x6_8_9 + y6_8_9 + z6_8_9 - 1.0 = 0;
subject to v7_8_9:
	-z7_8_8 - y7_7_9 - x6_8_9 + x7_8_9 + y7_8_9 + z7_8_9 - 1.0 = 0;
subject to v8_8_9:
	-z8_8_8 - y8_7_9 - x7_8_9 + x8_8_9 + y8_8_9 + z8_8_9 - 1.0 = 0;
subject to v9_8_9:
	-z9_8_8 - y9_7_9 - x8_8_9 + x9_8_9 + y9_8_9 + z9_8_9 - 1.0 = 0;
subject to v10_8_9:
	-x9_8_9 - z10_8_8 - y10_7_9 + y10_8_9 + z10_8_9 + y11_8_9 + z11_8_9 - 1.0 = 0;
subject to v1_9_9:
	-z1_9_8 - y1_8_9 + x1_9_9 + y1_9_9 + z1_9_9 + y0_9_9 + z0_9_9 - 1.0 = 0;
subject to v2_9_9:
	-z2_9_8 - y2_8_9 - x1_9_9 + x2_9_9 + y2_9_9 + z2_9_9 - 1.0 = 0;
subject to v3_9_9:
	-z3_9_8 - y3_8_9 - x2_9_9 + x3_9_9 + y3_9_9 + z3_9_9 - 1.0 = 0;
subject to v4_9_9:
	-z4_9_8 - y4_8_9 - x3_9_9 + x4_9_9 + y4_9_9 + z4_9_9 - 1.0 = 0;
subject to v5_9_9:
	-z5_9_8 - y5_8_9 - x4_9_9 + x5_9_9 + y5_9_9 + z5_9_9 - 1.0 = 0;
subject to v6_9_9:
	-z6_9_8 - y6_8_9 - x5_9_9 + x6_9_9 + y6_9_9 + z6_9_9 - 1.0 = 0;
subject to v7_9_9:
	-z7_9_8 - y7_8_9 - x6_9_9 + x7_9_9 + y7_9_9 + z7_9_9 - 1.0 = 0;
subject to v8_9_9:
	-z8_9_8 - y8_8_9 - x7_9_9 + x8_9_9 + y8_9_9 + z8_9_9 - 1.0 = 0;
subject to v9_9_9:
	-z9_9_8 - y9_8_9 - x8_9_9 + x9_9_9 + y9_9_9 + z9_9_9 - 1.0 = 0;
subject to v10_9_9:
	-x9_9_9 - z10_9_8 - y10_8_9 + y10_9_9 + z10_9_9 + y11_9_9 + z11_9_9 - 1.0 = 0;
subject to v1_10_9:
	-y1_9_9 - z1_10_8 + x1_10_9 + z1_10_9 + y0_10_9 + z0_10_9 + x1_11_9 + z1_11_9 - 
	1.0 = 0;
subject to v2_10_9:
	-y2_9_9 - z2_10_8 - x1_10_9 + x2_10_9 + z2_10_9 + x2_11_9 + z2_11_9 - 1.0 = 0;
subject to v3_10_9:
	-y3_9_9 - z3_10_8 - x2_10_9 + x3_10_9 + z3_10_9 + x3_11_9 + z3_11_9 - 1.0 = 0;
subject to v4_10_9:
	-y4_9_9 - z4_10_8 - x3_10_9 + x4_10_9 + z4_10_9 + x4_11_9 + z4_11_9 - 1.0 = 0;
subject to v5_10_9:
	-y5_9_9 - z5_10_8 - x4_10_9 + x5_10_9 + z5_10_9 + x5_11_9 + z5_11_9 - 1.0 = 0;
subject to v6_10_9:
	-y6_9_9 - z6_10_8 - x5_10_9 + x6_10_9 + z6_10_9 + x6_11_9 + z6_11_9 - 1.0 = 0;
subject to v7_10_9:
	-y7_9_9 - z7_10_8 - x6_10_9 + x7_10_9 + z7_10_9 + x7_11_9 + z7_11_9 - 1.0 = 0;
subject to v8_10_9:
	-y8_9_9 - z8_10_8 - x7_10_9 + x8_10_9 + z8_10_9 + x8_11_9 + z8_11_9 - 1.0 = 0;
subject to v9_10_9:
	-y9_9_9 - z9_10_8 - x8_10_9 + x9_10_9 + z9_10_9 + x9_11_9 + z9_11_9 - 1.0 = 0;
subject to v10_10_9:
	-y10_9_9 - x9_10_9 + y11_10_9 + z11_10_9 + x10_11_9 + z10_11_9 - 1.0 = 0;
subject to v1_1_10:
	-z1_1_9 + x1_1_10 + y1_1_10 + y0_1_10 + z0_1_10 + x1_0_10 + z1_0_10 + x1_1_11 + 
	y1_1_11 - 1.0 = 0;
subject to v2_1_10:
	-z2_1_9 - x1_1_10 + x2_1_10 + y2_1_10 + x2_0_10 + z2_0_10 + x2_1_11 + y2_1_11 - 
	1.0 = 0;
subject to v3_1_10:
	-z3_1_9 - x2_1_10 + x3_1_10 + y3_1_10 + x3_0_10 + z3_0_10 + x3_1_11 + y3_1_11 - 
	1.0 = 0;
subject to v4_1_10:
	-z4_1_9 - x3_1_10 + x4_1_10 + y4_1_10 + x4_0_10 + z4_0_10 + x4_1_11 + y4_1_11 - 
	1.0 = 0;
subject to v5_1_10:
	-z5_1_9 - x4_1_10 + x5_1_10 + y5_1_10 + x5_0_10 + z5_0_10 + x5_1_11 + y5_1_11 - 
	1.0 = 0;
subject to v6_1_10:
	-z6_1_9 - x5_1_10 + x6_1_10 + y6_1_10 + x6_0_10 + z6_0_10 + x6_1_11 + y6_1_11 - 
	1.0 = 0;
subject to v7_1_10:
	-z7_1_9 - x6_1_10 + x7_1_10 + y7_1_10 + x7_0_10 + z7_0_10 + x7_1_11 + y7_1_11 - 
	1.0 = 0;
subject to v8_1_10:
	-z8_1_9 - x7_1_10 + x8_1_10 + y8_1_10 + x8_0_10 + z8_0_10 + x8_1_11 + y8_1_11 - 
	1.0 = 0;
subject to v9_1_10:
	-z9_1_9 - x8_1_10 + x9_1_10 + y9_1_10 + x9_0_10 + z9_0_10 + x9_1_11 + y9_1_11 - 
	1.0 = 0;
subject to v10_1_10:
	-z10_1_9 - x9_1_10 + y11_1_10 + z11_1_10 + x10_0_10 + z10_0_10 + x10_1_11 + 
	y10_1_11 - 1.0 = 0;
subject to v1_2_10:
	-z1_2_9 - y1_1_10 + x1_2_10 + y1_2_10 + y0_2_10 + z0_2_10 + x1_2_11 + y1_2_11 - 
	1.0 = 0;
subject to v2_2_10:
	-z2_2_9 - y2_1_10 - x1_2_10 + x2_2_10 + y2_2_10 + x2_2_11 + y2_2_11 - 1.0 = 0;
subject to v3_2_10:
	-z3_2_9 - y3_1_10 - x2_2_10 + x3_2_10 + y3_2_10 + x3_2_11 + y3_2_11 - 1.0 = 0;
subject to v4_2_10:
	-z4_2_9 - y4_1_10 - x3_2_10 + x4_2_10 + y4_2_10 + x4_2_11 + y4_2_11 - 1.0 = 0;
subject to v5_2_10:
	-z5_2_9 - y5_1_10 - x4_2_10 + x5_2_10 + y5_2_10 + x5_2_11 + y5_2_11 - 1.0 = 0;
subject to v6_2_10:
	-z6_2_9 - y6_1_10 - x5_2_10 + x6_2_10 + y6_2_10 + x6_2_11 + y6_2_11 - 1.0 = 0;
subject to v7_2_10:
	-z7_2_9 - y7_1_10 - x6_2_10 + x7_2_10 + y7_2_10 + x7_2_11 + y7_2_11 - 1.0 = 0;
subject to v8_2_10:
	-z8_2_9 - y8_1_10 - x7_2_10 + x8_2_10 + y8_2_10 + x8_2_11 + y8_2_11 - 1.0 = 0;
subject to v9_2_10:
	-z9_2_9 - y9_1_10 - x8_2_10 + x9_2_10 + y9_2_10 + x9_2_11 + y9_2_11 - 1.0 = 0;
subject to v10_2_10:
	-z10_2_9 - x9_2_10 + y11_2_10 + z11_2_10 + x10_2_11 + y10_2_11 - 1.0 = 0;
subject to v1_3_10:
	-z1_3_9 - y1_2_10 + x1_3_10 + y1_3_10 + y0_3_10 + z0_3_10 + x1_3_11 + y1_3_11 - 
	1.0 = 0;
subject to v2_3_10:
	-z2_3_9 - y2_2_10 - x1_3_10 + x2_3_10 + y2_3_10 + x2_3_11 + y2_3_11 - 1.0 = 0;
subject to v3_3_10:
	-z3_3_9 - y3_2_10 - x2_3_10 + x3_3_10 + y3_3_10 + x3_3_11 + y3_3_11 - 1.0 = 0;
subject to v4_3_10:
	-z4_3_9 - y4_2_10 - x3_3_10 + x4_3_10 + y4_3_10 + x4_3_11 + y4_3_11 - 1.0 = 0;
subject to v5_3_10:
	-z5_3_9 - y5_2_10 - x4_3_10 + x5_3_10 + y5_3_10 + x5_3_11 + y5_3_11 - 1.0 = 0;
subject to v6_3_10:
	-z6_3_9 - y6_2_10 - x5_3_10 + x6_3_10 + y6_3_10 + x6_3_11 + y6_3_11 - 1.0 = 0;
subject to v7_3_10:
	-z7_3_9 - y7_2_10 - x6_3_10 + x7_3_10 + y7_3_10 + x7_3_11 + y7_3_11 - 1.0 = 0;
subject to v8_3_10:
	-z8_3_9 - y8_2_10 - x7_3_10 + x8_3_10 + y8_3_10 + x8_3_11 + y8_3_11 - 1.0 = 0;
subject to v9_3_10:
	-z9_3_9 - y9_2_10 - x8_3_10 + x9_3_10 + y9_3_10 + x9_3_11 + y9_3_11 - 1.0 = 0;
subject to v10_3_10:
	-z10_3_9 - x9_3_10 + y11_3_10 + z11_3_10 + x10_3_11 + y10_3_11 - 1.0 = 0;
subject to v1_4_10:
	-z1_4_9 - y1_3_10 + x1_4_10 + y1_4_10 + y0_4_10 + z0_4_10 + x1_4_11 + y1_4_11 - 
	1.0 = 0;
subject to v2_4_10:
	-z2_4_9 - y2_3_10 - x1_4_10 + x2_4_10 + y2_4_10 + x2_4_11 + y2_4_11 - 1.0 = 0;
subject to v3_4_10:
	-z3_4_9 - y3_3_10 - x2_4_10 + x3_4_10 + y3_4_10 + x3_4_11 + y3_4_11 - 1.0 = 0;
subject to v4_4_10:
	-z4_4_9 - y4_3_10 - x3_4_10 + x4_4_10 + y4_4_10 + x4_4_11 + y4_4_11 - 1.0 = 0;
subject to v5_4_10:
	-z5_4_9 - y5_3_10 - x4_4_10 + x5_4_10 + y5_4_10 + x5_4_11 + y5_4_11 - 1.0 = 0;
subject to v6_4_10:
	-z6_4_9 - y6_3_10 - x5_4_10 + x6_4_10 + y6_4_10 + x6_4_11 + y6_4_11 - 1.0 = 0;
subject to v7_4_10:
	-z7_4_9 - y7_3_10 - x6_4_10 + x7_4_10 + y7_4_10 + x7_4_11 + y7_4_11 - 1.0 = 0;
subject to v8_4_10:
	-z8_4_9 - y8_3_10 - x7_4_10 + x8_4_10 + y8_4_10 + x8_4_11 + y8_4_11 - 1.0 = 0;
subject to v9_4_10:
	-z9_4_9 - y9_3_10 - x8_4_10 + x9_4_10 + y9_4_10 + x9_4_11 + y9_4_11 - 1.0 = 0;
subject to v10_4_10:
	-z10_4_9 - x9_4_10 + y11_4_10 + z11_4_10 + x10_4_11 + y10_4_11 - 1.0 = 0;
subject to v1_5_10:
	-z1_5_9 - y1_4_10 + x1_5_10 + y1_5_10 + y0_5_10 + z0_5_10 + x1_5_11 + y1_5_11 - 
	1.0 = 0;
subject to v2_5_10:
	-z2_5_9 - y2_4_10 - x1_5_10 + x2_5_10 + y2_5_10 + x2_5_11 + y2_5_11 - 1.0 = 0;
subject to v3_5_10:
	-z3_5_9 - y3_4_10 - x2_5_10 + x3_5_10 + y3_5_10 + x3_5_11 + y3_5_11 - 1.0 = 0;
subject to v4_5_10:
	-z4_5_9 - y4_4_10 - x3_5_10 + x4_5_10 + y4_5_10 + x4_5_11 + y4_5_11 - 1.0 = 0;
subject to v5_5_10:
	-z5_5_9 - y5_4_10 - x4_5_10 + x5_5_10 + y5_5_10 + x5_5_11 + y5_5_11 - 1.0 = 0;
subject to v6_5_10:
	-z6_5_9 - y6_4_10 - x5_5_10 + x6_5_10 + y6_5_10 + x6_5_11 + y6_5_11 - 1.0 = 0;
subject to v7_5_10:
	-z7_5_9 - y7_4_10 - x6_5_10 + x7_5_10 + y7_5_10 + x7_5_11 + y7_5_11 - 1.0 = 0;
subject to v8_5_10:
	-z8_5_9 - y8_4_10 - x7_5_10 + x8_5_10 + y8_5_10 + x8_5_11 + y8_5_11 - 1.0 = 0;
subject to v9_5_10:
	-z9_5_9 - y9_4_10 - x8_5_10 + x9_5_10 + y9_5_10 + x9_5_11 + y9_5_11 - 1.0 = 0;
subject to v10_5_10:
	-z10_5_9 - x9_5_10 + y11_5_10 + z11_5_10 + x10_5_11 + y10_5_11 - 1.0 = 0;
subject to v1_6_10:
	-z1_6_9 - y1_5_10 + x1_6_10 + y1_6_10 + y0_6_10 + z0_6_10 + x1_6_11 + y1_6_11 - 
	1.0 = 0;
subject to v2_6_10:
	-z2_6_9 - y2_5_10 - x1_6_10 + x2_6_10 + y2_6_10 + x2_6_11 + y2_6_11 - 1.0 = 0;
subject to v3_6_10:
	-z3_6_9 - y3_5_10 - x2_6_10 + x3_6_10 + y3_6_10 + x3_6_11 + y3_6_11 - 1.0 = 0;
subject to v4_6_10:
	-z4_6_9 - y4_5_10 - x3_6_10 + x4_6_10 + y4_6_10 + x4_6_11 + y4_6_11 - 1.0 = 0;
subject to v5_6_10:
	-z5_6_9 - y5_5_10 - x4_6_10 + x5_6_10 + y5_6_10 + x5_6_11 + y5_6_11 - 1.0 = 0;
subject to v6_6_10:
	-z6_6_9 - y6_5_10 - x5_6_10 + x6_6_10 + y6_6_10 + x6_6_11 + y6_6_11 - 1.0 = 0;
subject to v7_6_10:
	-z7_6_9 - y7_5_10 - x6_6_10 + x7_6_10 + y7_6_10 + x7_6_11 + y7_6_11 - 1.0 = 0;
subject to v8_6_10:
	-z8_6_9 - y8_5_10 - x7_6_10 + x8_6_10 + y8_6_10 + x8_6_11 + y8_6_11 - 1.0 = 0;
subject to v9_6_10:
	-z9_6_9 - y9_5_10 - x8_6_10 + x9_6_10 + y9_6_10 + x9_6_11 + y9_6_11 - 1.0 = 0;
subject to v10_6_10:
	-z10_6_9 - x9_6_10 + y11_6_10 + z11_6_10 + x10_6_11 + y10_6_11 - 1.0 = 0;
subject to v1_7_10:
	-z1_7_9 - y1_6_10 + x1_7_10 + y1_7_10 + y0_7_10 + z0_7_10 + x1_7_11 + y1_7_11 - 
	1.0 = 0;
subject to v2_7_10:
	-z2_7_9 - y2_6_10 - x1_7_10 + x2_7_10 + y2_7_10 + x2_7_11 + y2_7_11 - 1.0 = 0;
subject to v3_7_10:
	-z3_7_9 - y3_6_10 - x2_7_10 + x3_7_10 + y3_7_10 + x3_7_11 + y3_7_11 - 1.0 = 0;
subject to v4_7_10:
	-z4_7_9 - y4_6_10 - x3_7_10 + x4_7_10 + y4_7_10 + x4_7_11 + y4_7_11 - 1.0 = 0;
subject to v5_7_10:
	-z5_7_9 - y5_6_10 - x4_7_10 + x5_7_10 + y5_7_10 + x5_7_11 + y5_7_11 - 1.0 = 0;
subject to v6_7_10:
	-z6_7_9 - y6_6_10 - x5_7_10 + x6_7_10 + y6_7_10 + x6_7_11 + y6_7_11 - 1.0 = 0;
subject to v7_7_10:
	-z7_7_9 - y7_6_10 - x6_7_10 + x7_7_10 + y7_7_10 + x7_7_11 + y7_7_11 - 1.0 = 0;
subject to v8_7_10:
	-z8_7_9 - y8_6_10 - x7_7_10 + x8_7_10 + y8_7_10 + x8_7_11 + y8_7_11 - 1.0 = 0;
subject to v9_7_10:
	-z9_7_9 - y9_6_10 - x8_7_10 + x9_7_10 + y9_7_10 + x9_7_11 + y9_7_11 - 1.0 = 0;
subject to v10_7_10:
	-z10_7_9 - x9_7_10 + y11_7_10 + z11_7_10 + x10_7_11 + y10_7_11 - 1.0 = 0;
subject to v1_8_10:
	-z1_8_9 - y1_7_10 + x1_8_10 + y1_8_10 + y0_8_10 + z0_8_10 + x1_8_11 + y1_8_11 - 
	1.0 = 0;
subject to v2_8_10:
	-z2_8_9 - y2_7_10 - x1_8_10 + x2_8_10 + y2_8_10 + x2_8_11 + y2_8_11 - 1.0 = 0;
subject to v3_8_10:
	-z3_8_9 - y3_7_10 - x2_8_10 + x3_8_10 + y3_8_10 + x3_8_11 + y3_8_11 - 1.0 = 0;
subject to v4_8_10:
	-z4_8_9 - y4_7_10 - x3_8_10 + x4_8_10 + y4_8_10 + x4_8_11 + y4_8_11 - 1.0 = 0;
subject to v5_8_10:
	-z5_8_9 - y5_7_10 - x4_8_10 + x5_8_10 + y5_8_10 + x5_8_11 + y5_8_11 - 1.0 = 0;
subject to v6_8_10:
	-z6_8_9 - y6_7_10 - x5_8_10 + x6_8_10 + y6_8_10 + x6_8_11 + y6_8_11 - 1.0 = 0;
subject to v7_8_10:
	-z7_8_9 - y7_7_10 - x6_8_10 + x7_8_10 + y7_8_10 + x7_8_11 + y7_8_11 - 1.0 = 0;
subject to v8_8_10:
	-z8_8_9 - y8_7_10 - x7_8_10 + x8_8_10 + y8_8_10 + x8_8_11 + y8_8_11 - 1.0 = 0;
subject to v9_8_10:
	-z9_8_9 - y9_7_10 - x8_8_10 + x9_8_10 + y9_8_10 + x9_8_11 + y9_8_11 - 1.0 = 0;
subject to v10_8_10:
	-z10_8_9 - x9_8_10 + y11_8_10 + z11_8_10 + x10_8_11 + y10_8_11 - 1.0 = 0;
subject to v1_9_10:
	-z1_9_9 - y1_8_10 + x1_9_10 + y1_9_10 + y0_9_10 + z0_9_10 + x1_9_11 + y1_9_11 - 
	1.0 = 0;
subject to v2_9_10:
	-z2_9_9 - y2_8_10 - x1_9_10 + x2_9_10 + y2_9_10 + x2_9_11 + y2_9_11 - 1.0 = 0;
subject to v3_9_10:
	-z3_9_9 - y3_8_10 - x2_9_10 + x3_9_10 + y3_9_10 + x3_9_11 + y3_9_11 - 1.0 = 0;
subject to v4_9_10:
	-z4_9_9 - y4_8_10 - x3_9_10 + x4_9_10 + y4_9_10 + x4_9_11 + y4_9_11 - 1.0 = 0;
subject to v5_9_10:
	-z5_9_9 - y5_8_10 - x4_9_10 + x5_9_10 + y5_9_10 + x5_9_11 + y5_9_11 - 1.0 = 0;
subject to v6_9_10:
	-z6_9_9 - y6_8_10 - x5_9_10 + x6_9_10 + y6_9_10 + x6_9_11 + y6_9_11 - 1.0 = 0;
subject to v7_9_10:
	-z7_9_9 - y7_8_10 - x6_9_10 + x7_9_10 + y7_9_10 + x7_9_11 + y7_9_11 - 1.0 = 0;
subject to v8_9_10:
	-z8_9_9 - y8_8_10 - x7_9_10 + x8_9_10 + y8_9_10 + x8_9_11 + y8_9_11 - 1.0 = 0;
subject to v9_9_10:
	-z9_9_9 - y9_8_10 - x8_9_10 + x9_9_10 + y9_9_10 + x9_9_11 + y9_9_11 - 1.0 = 0;
subject to v10_9_10:
	-z10_9_9 - x9_9_10 + y11_9_10 + z11_9_10 + x10_9_11 + y10_9_11 - 1.0 = 0;
subject to v1_10_10:
	-z1_10_9 - y1_9_10 + y0_10_10 + z0_10_10 + x1_11_10 + z1_11_10 + x1_10_11 + 
	y1_10_11 - 1.0 = 0;
subject to v2_10_10:
	-z2_10_9 - y2_9_10 + x2_11_10 + z2_11_10 + x2_10_11 + y2_10_11 - 1.0 = 0;
subject to v3_10_10:
	-z3_10_9 - y3_9_10 + x3_11_10 + z3_11_10 + x3_10_11 + y3_10_11 - 1.0 = 0;
subject to v4_10_10:
	-z4_10_9 - y4_9_10 + x4_11_10 + z4_11_10 + x4_10_11 + y4_10_11 - 1.0 = 0;
subject to v5_10_10:
	-z5_10_9 - y5_9_10 + x5_11_10 + z5_11_10 + x5_10_11 + y5_10_11 - 1.0 = 0;
subject to v6_10_10:
	-z6_10_9 - y6_9_10 + x6_11_10 + z6_11_10 + x6_10_11 + y6_10_11 - 1.0 = 0;
subject to v7_10_10:
	-z7_10_9 - y7_9_10 + x7_11_10 + z7_11_10 + x7_10_11 + y7_10_11 - 1.0 = 0;
subject to v8_10_10:
	-z8_10_9 - y8_9_10 + x8_11_10 + z8_11_10 + x8_10_11 + y8_10_11 - 1.0 = 0;
subject to v9_10_10:
	-z9_10_9 - y9_9_10 + x9_11_10 + z9_11_10 + x9_10_11 + y9_10_11 - 1.0 = 0;
subject to v10_10_10:
	y11_10_10 + z11_10_10 + x10_11_10 + z10_11_10 + x10_10_11 + y10_10_11 - 1.0 = 0;

solve;
	display x1_1_1;
	display y1_1_1;
	display z1_1_1;
	display x2_1_1;
	display y2_1_1;
	display z2_1_1;
	display x3_1_1;
	display y3_1_1;
	display z3_1_1;
	display x4_1_1;
	display y4_1_1;
	display z4_1_1;
	display x5_1_1;
	display y5_1_1;
	display z5_1_1;
	display x6_1_1;
	display y6_1_1;
	display z6_1_1;
	display x7_1_1;
	display y7_1_1;
	display z7_1_1;
	display x8_1_1;
	display y8_1_1;
	display z8_1_1;
	display x9_1_1;
	display y9_1_1;
	display z9_1_1;
	display x1_2_1;
	display y1_2_1;
	display z1_2_1;
	display x2_2_1;
	display y2_2_1;
	display z2_2_1;
	display x3_2_1;
	display y3_2_1;
	display z3_2_1;
	display x4_2_1;
	display y4_2_1;
	display z4_2_1;
	display x5_2_1;
	display y5_2_1;
	display z5_2_1;
	display x6_2_1;
	display y6_2_1;
	display z6_2_1;
	display x7_2_1;
	display y7_2_1;
	display z7_2_1;
	display x8_2_1;
	display y8_2_1;
	display z8_2_1;
	display x9_2_1;
	display y9_2_1;
	display z9_2_1;
	display x1_3_1;
	display y1_3_1;
	display z1_3_1;
	display x2_3_1;
	display y2_3_1;
	display z2_3_1;
	display x3_3_1;
	display y3_3_1;
	display z3_3_1;
	display x4_3_1;
	display y4_3_1;
	display z4_3_1;
	display x5_3_1;
	display y5_3_1;
	display z5_3_1;
	display x6_3_1;
	display y6_3_1;
	display z6_3_1;
	display x7_3_1;
	display y7_3_1;
	display z7_3_1;
	display x8_3_1;
	display y8_3_1;
	display z8_3_1;
	display x9_3_1;
	display y9_3_1;
	display z9_3_1;
	display x1_4_1;
	display y1_4_1;
	display z1_4_1;
	display x2_4_1;
	display y2_4_1;
	display z2_4_1;
	display x3_4_1;
	display y3_4_1;
	display z3_4_1;
	display x4_4_1;
	display y4_4_1;
	display z4_4_1;
	display x5_4_1;
	display y5_4_1;
	display z5_4_1;
	display x6_4_1;
	display y6_4_1;
	display z6_4_1;
	display x7_4_1;
	display y7_4_1;
	display z7_4_1;
	display x8_4_1;
	display y8_4_1;
	display z8_4_1;
	display x9_4_1;
	display y9_4_1;
	display z9_4_1;
	display x1_5_1;
	display y1_5_1;
	display z1_5_1;
	display x2_5_1;
	display y2_5_1;
	display z2_5_1;
	display x3_5_1;
	display y3_5_1;
	display z3_5_1;
	display x4_5_1;
	display y4_5_1;
	display z4_5_1;
	display x5_5_1;
	display y5_5_1;
	display z5_5_1;
	display x6_5_1;
	display y6_5_1;
	display z6_5_1;
	display x7_5_1;
	display y7_5_1;
	display z7_5_1;
	display x8_5_1;
	display y8_5_1;
	display z8_5_1;
	display x9_5_1;
	display y9_5_1;
	display z9_5_1;
	display x1_6_1;
	display y1_6_1;
	display z1_6_1;
	display x2_6_1;
	display y2_6_1;
	display z2_6_1;
	display x3_6_1;
	display y3_6_1;
	display z3_6_1;
	display x4_6_1;
	display y4_6_1;
	display z4_6_1;
	display x5_6_1;
	display y5_6_1;
	display z5_6_1;
	display x6_6_1;
	display y6_6_1;
	display z6_6_1;
	display x7_6_1;
	display y7_6_1;
	display z7_6_1;
	display x8_6_1;
	display y8_6_1;
	display z8_6_1;
	display x9_6_1;
	display y9_6_1;
	display z9_6_1;
	display x1_7_1;
	display y1_7_1;
	display z1_7_1;
	display x2_7_1;
	display y2_7_1;
	display z2_7_1;
	display x3_7_1;
	display y3_7_1;
	display z3_7_1;
	display x4_7_1;
	display y4_7_1;
	display z4_7_1;
	display x5_7_1;
	display y5_7_1;
	display z5_7_1;
	display x6_7_1;
	display y6_7_1;
	display z6_7_1;
	display x7_7_1;
	display y7_7_1;
	display z7_7_1;
	display x8_7_1;
	display y8_7_1;
	display z8_7_1;
	display x9_7_1;
	display y9_7_1;
	display z9_7_1;
	display x1_8_1;
	display y1_8_1;
	display z1_8_1;
	display x2_8_1;
	display y2_8_1;
	display z2_8_1;
	display x3_8_1;
	display y3_8_1;
	display z3_8_1;
	display x4_8_1;
	display y4_8_1;
	display z4_8_1;
	display x5_8_1;
	display y5_8_1;
	display z5_8_1;
	display x6_8_1;
	display y6_8_1;
	display z6_8_1;
	display x7_8_1;
	display y7_8_1;
	display z7_8_1;
	display x8_8_1;
	display y8_8_1;
	display z8_8_1;
	display x9_8_1;
	display y9_8_1;
	display z9_8_1;
	display x1_9_1;
	display y1_9_1;
	display z1_9_1;
	display x2_9_1;
	display y2_9_1;
	display z2_9_1;
	display x3_9_1;
	display y3_9_1;
	display z3_9_1;
	display x4_9_1;
	display y4_9_1;
	display z4_9_1;
	display x5_9_1;
	display y5_9_1;
	display z5_9_1;
	display x6_9_1;
	display y6_9_1;
	display z6_9_1;
	display x7_9_1;
	display y7_9_1;
	display z7_9_1;
	display x8_9_1;
	display y8_9_1;
	display z8_9_1;
	display x9_9_1;
	display y9_9_1;
	display z9_9_1;
	display x1_1_2;
	display y1_1_2;
	display z1_1_2;
	display x2_1_2;
	display y2_1_2;
	display z2_1_2;
	display x3_1_2;
	display y3_1_2;
	display z3_1_2;
	display x4_1_2;
	display y4_1_2;
	display z4_1_2;
	display x5_1_2;
	display y5_1_2;
	display z5_1_2;
	display x6_1_2;
	display y6_1_2;
	display z6_1_2;
	display x7_1_2;
	display y7_1_2;
	display z7_1_2;
	display x8_1_2;
	display y8_1_2;
	display z8_1_2;
	display x9_1_2;
	display y9_1_2;
	display z9_1_2;
	display x1_2_2;
	display y1_2_2;
	display z1_2_2;
	display x2_2_2;
	display y2_2_2;
	display z2_2_2;
	display x3_2_2;
	display y3_2_2;
	display z3_2_2;
	display x4_2_2;
	display y4_2_2;
	display z4_2_2;
	display x5_2_2;
	display y5_2_2;
	display z5_2_2;
	display x6_2_2;
	display y6_2_2;
	display z6_2_2;
	display x7_2_2;
	display y7_2_2;
	display z7_2_2;
	display x8_2_2;
	display y8_2_2;
	display z8_2_2;
	display x9_2_2;
	display y9_2_2;
	display z9_2_2;
	display x1_3_2;
	display y1_3_2;
	display z1_3_2;
	display x2_3_2;
	display y2_3_2;
	display z2_3_2;
	display x3_3_2;
	display y3_3_2;
	display z3_3_2;
	display x4_3_2;
	display y4_3_2;
	display z4_3_2;
	display x5_3_2;
	display y5_3_2;
	display z5_3_2;
	display x6_3_2;
	display y6_3_2;
	display z6_3_2;
	display x7_3_2;
	display y7_3_2;
	display z7_3_2;
	display x8_3_2;
	display y8_3_2;
	display z8_3_2;
	display x9_3_2;
	display y9_3_2;
	display z9_3_2;
	display x1_4_2;
	display y1_4_2;
	display z1_4_2;
	display x2_4_2;
	display y2_4_2;
	display z2_4_2;
	display x3_4_2;
	display y3_4_2;
	display z3_4_2;
	display x4_4_2;
	display y4_4_2;
	display z4_4_2;
	display x5_4_2;
	display y5_4_2;
	display z5_4_2;
	display x6_4_2;
	display y6_4_2;
	display z6_4_2;
	display x7_4_2;
	display y7_4_2;
	display z7_4_2;
	display x8_4_2;
	display y8_4_2;
	display z8_4_2;
	display x9_4_2;
	display y9_4_2;
	display z9_4_2;
	display x1_5_2;
	display y1_5_2;
	display z1_5_2;
	display x2_5_2;
	display y2_5_2;
	display z2_5_2;
	display x3_5_2;
	display y3_5_2;
	display z3_5_2;
	display x4_5_2;
	display y4_5_2;
	display z4_5_2;
	display x5_5_2;
	display y5_5_2;
	display z5_5_2;
	display x6_5_2;
	display y6_5_2;
	display z6_5_2;
	display x7_5_2;
	display y7_5_2;
	display z7_5_2;
	display x8_5_2;
	display y8_5_2;
	display z8_5_2;
	display x9_5_2;
	display y9_5_2;
	display z9_5_2;
	display x1_6_2;
	display y1_6_2;
	display z1_6_2;
	display x2_6_2;
	display y2_6_2;
	display z2_6_2;
	display x3_6_2;
	display y3_6_2;
	display z3_6_2;
	display x4_6_2;
	display y4_6_2;
	display z4_6_2;
	display x5_6_2;
	display y5_6_2;
	display z5_6_2;
	display x6_6_2;
	display y6_6_2;
	display z6_6_2;
	display x7_6_2;
	display y7_6_2;
	display z7_6_2;
	display x8_6_2;
	display y8_6_2;
	display z8_6_2;
	display x9_6_2;
	display y9_6_2;
	display z9_6_2;
	display x1_7_2;
	display y1_7_2;
	display z1_7_2;
	display x2_7_2;
	display y2_7_2;
	display z2_7_2;
	display x3_7_2;
	display y3_7_2;
	display z3_7_2;
	display x4_7_2;
	display y4_7_2;
	display z4_7_2;
	display x5_7_2;
	display y5_7_2;
	display z5_7_2;
	display x6_7_2;
	display y6_7_2;
	display z6_7_2;
	display x7_7_2;
	display y7_7_2;
	display z7_7_2;
	display x8_7_2;
	display y8_7_2;
	display z8_7_2;
	display x9_7_2;
	display y9_7_2;
	display z9_7_2;
	display x1_8_2;
	display y1_8_2;
	display z1_8_2;
	display x2_8_2;
	display y2_8_2;
	display z2_8_2;
	display x3_8_2;
	display y3_8_2;
	display z3_8_2;
	display x4_8_2;
	display y4_8_2;
	display z4_8_2;
	display x5_8_2;
	display y5_8_2;
	display z5_8_2;
	display x6_8_2;
	display y6_8_2;
	display z6_8_2;
	display x7_8_2;
	display y7_8_2;
	display z7_8_2;
	display x8_8_2;
	display y8_8_2;
	display z8_8_2;
	display x9_8_2;
	display y9_8_2;
	display z9_8_2;
	display x1_9_2;
	display y1_9_2;
	display z1_9_2;
	display x2_9_2;
	display y2_9_2;
	display z2_9_2;
	display x3_9_2;
	display y3_9_2;
	display z3_9_2;
	display x4_9_2;
	display y4_9_2;
	display z4_9_2;
	display x5_9_2;
	display y5_9_2;
	display z5_9_2;
	display x6_9_2;
	display y6_9_2;
	display z6_9_2;
	display x7_9_2;
	display y7_9_2;
	display z7_9_2;
	display x8_9_2;
	display y8_9_2;
	display z8_9_2;
	display x9_9_2;
	display y9_9_2;
	display z9_9_2;
	display x1_1_3;
	display y1_1_3;
	display z1_1_3;
	display x2_1_3;
	display y2_1_3;
	display z2_1_3;
	display x3_1_3;
	display y3_1_3;
	display z3_1_3;
	display x4_1_3;
	display y4_1_3;
	display z4_1_3;
	display x5_1_3;
	display y5_1_3;
	display z5_1_3;
	display x6_1_3;
	display y6_1_3;
	display z6_1_3;
	display x7_1_3;
	display y7_1_3;
	display z7_1_3;
	display x8_1_3;
	display y8_1_3;
	display z8_1_3;
	display x9_1_3;
	display y9_1_3;
	display z9_1_3;
	display x1_2_3;
	display y1_2_3;
	display z1_2_3;
	display x2_2_3;
	display y2_2_3;
	display z2_2_3;
	display x3_2_3;
	display y3_2_3;
	display z3_2_3;
	display x4_2_3;
	display y4_2_3;
	display z4_2_3;
	display x5_2_3;
	display y5_2_3;
	display z5_2_3;
	display x6_2_3;
	display y6_2_3;
	display z6_2_3;
	display x7_2_3;
	display y7_2_3;
	display z7_2_3;
	display x8_2_3;
	display y8_2_3;
	display z8_2_3;
	display x9_2_3;
	display y9_2_3;
	display z9_2_3;
	display x1_3_3;
	display y1_3_3;
	display z1_3_3;
	display x2_3_3;
	display y2_3_3;
	display z2_3_3;
	display x3_3_3;
	display y3_3_3;
	display z3_3_3;
	display x4_3_3;
	display y4_3_3;
	display z4_3_3;
	display x5_3_3;
	display y5_3_3;
	display z5_3_3;
	display x6_3_3;
	display y6_3_3;
	display z6_3_3;
	display x7_3_3;
	display y7_3_3;
	display z7_3_3;
	display x8_3_3;
	display y8_3_3;
	display z8_3_3;
	display x9_3_3;
	display y9_3_3;
	display z9_3_3;
	display x1_4_3;
	display y1_4_3;
	display z1_4_3;
	display x2_4_3;
	display y2_4_3;
	display z2_4_3;
	display x3_4_3;
	display y3_4_3;
	display z3_4_3;
	display x4_4_3;
	display y4_4_3;
	display z4_4_3;
	display x5_4_3;
	display y5_4_3;
	display z5_4_3;
	display x6_4_3;
	display y6_4_3;
	display z6_4_3;
	display x7_4_3;
	display y7_4_3;
	display z7_4_3;
	display x8_4_3;
	display y8_4_3;
	display z8_4_3;
	display x9_4_3;
	display y9_4_3;
	display z9_4_3;
	display x1_5_3;
	display y1_5_3;
	display z1_5_3;
	display x2_5_3;
	display y2_5_3;
	display z2_5_3;
	display x3_5_3;
	display y3_5_3;
	display z3_5_3;
	display x4_5_3;
	display y4_5_3;
	display z4_5_3;
	display x5_5_3;
	display y5_5_3;
	display z5_5_3;
	display x6_5_3;
	display y6_5_3;
	display z6_5_3;
	display x7_5_3;
	display y7_5_3;
	display z7_5_3;
	display x8_5_3;
	display y8_5_3;
	display z8_5_3;
	display x9_5_3;
	display y9_5_3;
	display z9_5_3;
	display x1_6_3;
	display y1_6_3;
	display z1_6_3;
	display x2_6_3;
	display y2_6_3;
	display z2_6_3;
	display x3_6_3;
	display y3_6_3;
	display z3_6_3;
	display x4_6_3;
	display y4_6_3;
	display z4_6_3;
	display x5_6_3;
	display y5_6_3;
	display z5_6_3;
	display x6_6_3;
	display y6_6_3;
	display z6_6_3;
	display x7_6_3;
	display y7_6_3;
	display z7_6_3;
	display x8_6_3;
	display y8_6_3;
	display z8_6_3;
	display x9_6_3;
	display y9_6_3;
	display z9_6_3;
	display x1_7_3;
	display y1_7_3;
	display z1_7_3;
	display x2_7_3;
	display y2_7_3;
	display z2_7_3;
	display x3_7_3;
	display y3_7_3;
	display z3_7_3;
	display x4_7_3;
	display y4_7_3;
	display z4_7_3;
	display x5_7_3;
	display y5_7_3;
	display z5_7_3;
	display x6_7_3;
	display y6_7_3;
	display z6_7_3;
	display x7_7_3;
	display y7_7_3;
	display z7_7_3;
	display x8_7_3;
	display y8_7_3;
	display z8_7_3;
	display x9_7_3;
	display y9_7_3;
	display z9_7_3;
	display x1_8_3;
	display y1_8_3;
	display z1_8_3;
	display x2_8_3;
	display y2_8_3;
	display z2_8_3;
	display x3_8_3;
	display y3_8_3;
	display z3_8_3;
	display x4_8_3;
	display y4_8_3;
	display z4_8_3;
	display x5_8_3;
	display y5_8_3;
	display z5_8_3;
	display x6_8_3;
	display y6_8_3;
	display z6_8_3;
	display x7_8_3;
	display y7_8_3;
	display z7_8_3;
	display x8_8_3;
	display y8_8_3;
	display z8_8_3;
	display x9_8_3;
	display y9_8_3;
	display z9_8_3;
	display x1_9_3;
	display y1_9_3;
	display z1_9_3;
	display x2_9_3;
	display y2_9_3;
	display z2_9_3;
	display x3_9_3;
	display y3_9_3;
	display z3_9_3;
	display x4_9_3;
	display y4_9_3;
	display z4_9_3;
	display x5_9_3;
	display y5_9_3;
	display z5_9_3;
	display x6_9_3;
	display y6_9_3;
	display z6_9_3;
	display x7_9_3;
	display y7_9_3;
	display z7_9_3;
	display x8_9_3;
	display y8_9_3;
	display z8_9_3;
	display x9_9_3;
	display y9_9_3;
	display z9_9_3;
	display x1_1_4;
	display y1_1_4;
	display z1_1_4;
	display x2_1_4;
	display y2_1_4;
	display z2_1_4;
	display x3_1_4;
	display y3_1_4;
	display z3_1_4;
	display x4_1_4;
	display y4_1_4;
	display z4_1_4;
	display x5_1_4;
	display y5_1_4;
	display z5_1_4;
	display x6_1_4;
	display y6_1_4;
	display z6_1_4;
	display x7_1_4;
	display y7_1_4;
	display z7_1_4;
	display x8_1_4;
	display y8_1_4;
	display z8_1_4;
	display x9_1_4;
	display y9_1_4;
	display z9_1_4;
	display x1_2_4;
	display y1_2_4;
	display z1_2_4;
	display x2_2_4;
	display y2_2_4;
	display z2_2_4;
	display x3_2_4;
	display y3_2_4;
	display z3_2_4;
	display x4_2_4;
	display y4_2_4;
	display z4_2_4;
	display x5_2_4;
	display y5_2_4;
	display z5_2_4;
	display x6_2_4;
	display y6_2_4;
	display z6_2_4;
	display x7_2_4;
	display y7_2_4;
	display z7_2_4;
	display x8_2_4;
	display y8_2_4;
	display z8_2_4;
	display x9_2_4;
	display y9_2_4;
	display z9_2_4;
	display x1_3_4;
	display y1_3_4;
	display z1_3_4;
	display x2_3_4;
	display y2_3_4;
	display z2_3_4;
	display x3_3_4;
	display y3_3_4;
	display z3_3_4;
	display x4_3_4;
	display y4_3_4;
	display z4_3_4;
	display x5_3_4;
	display y5_3_4;
	display z5_3_4;
	display x6_3_4;
	display y6_3_4;
	display z6_3_4;
	display x7_3_4;
	display y7_3_4;
	display z7_3_4;
	display x8_3_4;
	display y8_3_4;
	display z8_3_4;
	display x9_3_4;
	display y9_3_4;
	display z9_3_4;
	display x1_4_4;
	display y1_4_4;
	display z1_4_4;
	display x2_4_4;
	display y2_4_4;
	display z2_4_4;
	display x3_4_4;
	display y3_4_4;
	display z3_4_4;
	display x4_4_4;
	display y4_4_4;
	display z4_4_4;
	display x5_4_4;
	display y5_4_4;
	display z5_4_4;
	display x6_4_4;
	display y6_4_4;
	display z6_4_4;
	display x7_4_4;
	display y7_4_4;
	display z7_4_4;
	display x8_4_4;
	display y8_4_4;
	display z8_4_4;
	display x9_4_4;
	display y9_4_4;
	display z9_4_4;
	display x1_5_4;
	display y1_5_4;
	display z1_5_4;
	display x2_5_4;
	display y2_5_4;
	display z2_5_4;
	display x3_5_4;
	display y3_5_4;
	display z3_5_4;
	display x4_5_4;
	display y4_5_4;
	display z4_5_4;
	display x5_5_4;
	display y5_5_4;
	display z5_5_4;
	display x6_5_4;
	display y6_5_4;
	display z6_5_4;
	display x7_5_4;
	display y7_5_4;
	display z7_5_4;
	display x8_5_4;
	display y8_5_4;
	display z8_5_4;
	display x9_5_4;
	display y9_5_4;
	display z9_5_4;
	display x1_6_4;
	display y1_6_4;
	display z1_6_4;
	display x2_6_4;
	display y2_6_4;
	display z2_6_4;
	display x3_6_4;
	display y3_6_4;
	display z3_6_4;
	display x4_6_4;
	display y4_6_4;
	display z4_6_4;
	display x5_6_4;
	display y5_6_4;
	display z5_6_4;
	display x6_6_4;
	display y6_6_4;
	display z6_6_4;
	display x7_6_4;
	display y7_6_4;
	display z7_6_4;
	display x8_6_4;
	display y8_6_4;
	display z8_6_4;
	display x9_6_4;
	display y9_6_4;
	display z9_6_4;
	display x1_7_4;
	display y1_7_4;
	display z1_7_4;
	display x2_7_4;
	display y2_7_4;
	display z2_7_4;
	display x3_7_4;
	display y3_7_4;
	display z3_7_4;
	display x4_7_4;
	display y4_7_4;
	display z4_7_4;
	display x5_7_4;
	display y5_7_4;
	display z5_7_4;
	display x6_7_4;
	display y6_7_4;
	display z6_7_4;
	display x7_7_4;
	display y7_7_4;
	display z7_7_4;
	display x8_7_4;
	display y8_7_4;
	display z8_7_4;
	display x9_7_4;
	display y9_7_4;
	display z9_7_4;
	display x1_8_4;
	display y1_8_4;
	display z1_8_4;
	display x2_8_4;
	display y2_8_4;
	display z2_8_4;
	display x3_8_4;
	display y3_8_4;
	display z3_8_4;
	display x4_8_4;
	display y4_8_4;
	display z4_8_4;
	display x5_8_4;
	display y5_8_4;
	display z5_8_4;
	display x6_8_4;
	display y6_8_4;
	display z6_8_4;
	display x7_8_4;
	display y7_8_4;
	display z7_8_4;
	display x8_8_4;
	display y8_8_4;
	display z8_8_4;
	display x9_8_4;
	display y9_8_4;
	display z9_8_4;
	display x1_9_4;
	display y1_9_4;
	display z1_9_4;
	display x2_9_4;
	display y2_9_4;
	display z2_9_4;
	display x3_9_4;
	display y3_9_4;
	display z3_9_4;
	display x4_9_4;
	display y4_9_4;
	display z4_9_4;
	display x5_9_4;
	display y5_9_4;
	display z5_9_4;
	display x6_9_4;
	display y6_9_4;
	display z6_9_4;
	display x7_9_4;
	display y7_9_4;
	display z7_9_4;
	display x8_9_4;
	display y8_9_4;
	display z8_9_4;
	display x9_9_4;
	display y9_9_4;
	display z9_9_4;
	display x1_1_5;
	display y1_1_5;
	display z1_1_5;
	display x2_1_5;
	display y2_1_5;
	display z2_1_5;
	display x3_1_5;
	display y3_1_5;
	display z3_1_5;
	display x4_1_5;
	display y4_1_5;
	display z4_1_5;
	display x5_1_5;
	display y5_1_5;
	display z5_1_5;
	display x6_1_5;
	display y6_1_5;
	display z6_1_5;
	display x7_1_5;
	display y7_1_5;
	display z7_1_5;
	display x8_1_5;
	display y8_1_5;
	display z8_1_5;
	display x9_1_5;
	display y9_1_5;
	display z9_1_5;
	display x1_2_5;
	display y1_2_5;
	display z1_2_5;
	display x2_2_5;
	display y2_2_5;
	display z2_2_5;
	display x3_2_5;
	display y3_2_5;
	display z3_2_5;
	display x4_2_5;
	display y4_2_5;
	display z4_2_5;
	display x5_2_5;
	display y5_2_5;
	display z5_2_5;
	display x6_2_5;
	display y6_2_5;
	display z6_2_5;
	display x7_2_5;
	display y7_2_5;
	display z7_2_5;
	display x8_2_5;
	display y8_2_5;
	display z8_2_5;
	display x9_2_5;
	display y9_2_5;
	display z9_2_5;
	display x1_3_5;
	display y1_3_5;
	display z1_3_5;
	display x2_3_5;
	display y2_3_5;
	display z2_3_5;
	display x3_3_5;
	display y3_3_5;
	display z3_3_5;
	display x4_3_5;
	display y4_3_5;
	display z4_3_5;
	display x5_3_5;
	display y5_3_5;
	display z5_3_5;
	display x6_3_5;
	display y6_3_5;
	display z6_3_5;
	display x7_3_5;
	display y7_3_5;
	display z7_3_5;
	display x8_3_5;
	display y8_3_5;
	display z8_3_5;
	display x9_3_5;
	display y9_3_5;
	display z9_3_5;
	display x1_4_5;
	display y1_4_5;
	display z1_4_5;
	display x2_4_5;
	display y2_4_5;
	display z2_4_5;
	display x3_4_5;
	display y3_4_5;
	display z3_4_5;
	display x4_4_5;
	display y4_4_5;
	display z4_4_5;
	display x5_4_5;
	display y5_4_5;
	display z5_4_5;
	display x6_4_5;
	display y6_4_5;
	display z6_4_5;
	display x7_4_5;
	display y7_4_5;
	display z7_4_5;
	display x8_4_5;
	display y8_4_5;
	display z8_4_5;
	display x9_4_5;
	display y9_4_5;
	display z9_4_5;
	display x1_5_5;
	display y1_5_5;
	display z1_5_5;
	display x2_5_5;
	display y2_5_5;
	display z2_5_5;
	display x3_5_5;
	display y3_5_5;
	display z3_5_5;
	display x4_5_5;
	display y4_5_5;
	display z4_5_5;
	display x5_5_5;
	display y5_5_5;
	display z5_5_5;
	display x6_5_5;
	display y6_5_5;
	display z6_5_5;
	display x7_5_5;
	display y7_5_5;
	display z7_5_5;
	display x8_5_5;
	display y8_5_5;
	display z8_5_5;
	display x9_5_5;
	display y9_5_5;
	display z9_5_5;
	display x1_6_5;
	display y1_6_5;
	display z1_6_5;
	display x2_6_5;
	display y2_6_5;
	display z2_6_5;
	display x3_6_5;
	display y3_6_5;
	display z3_6_5;
	display x4_6_5;
	display y4_6_5;
	display z4_6_5;
	display x5_6_5;
	display y5_6_5;
	display z5_6_5;
	display x6_6_5;
	display y6_6_5;
	display z6_6_5;
	display x7_6_5;
	display y7_6_5;
	display z7_6_5;
	display x8_6_5;
	display y8_6_5;
	display z8_6_5;
	display x9_6_5;
	display y9_6_5;
	display z9_6_5;
	display x1_7_5;
	display y1_7_5;
	display z1_7_5;
	display x2_7_5;
	display y2_7_5;
	display z2_7_5;
	display x3_7_5;
	display y3_7_5;
	display z3_7_5;
	display x4_7_5;
	display y4_7_5;
	display z4_7_5;
	display x5_7_5;
	display y5_7_5;
	display z5_7_5;
	display x6_7_5;
	display y6_7_5;
	display z6_7_5;
	display x7_7_5;
	display y7_7_5;
	display z7_7_5;
	display x8_7_5;
	display y8_7_5;
	display z8_7_5;
	display x9_7_5;
	display y9_7_5;
	display z9_7_5;
	display x1_8_5;
	display y1_8_5;
	display z1_8_5;
	display x2_8_5;
	display y2_8_5;
	display z2_8_5;
	display x3_8_5;
	display y3_8_5;
	display z3_8_5;
	display x4_8_5;
	display y4_8_5;
	display z4_8_5;
	display x5_8_5;
	display y5_8_5;
	display z5_8_5;
	display x6_8_5;
	display y6_8_5;
	display z6_8_5;
	display x7_8_5;
	display y7_8_5;
	display z7_8_5;
	display x8_8_5;
	display y8_8_5;
	display z8_8_5;
	display x9_8_5;
	display y9_8_5;
	display z9_8_5;
	display x1_9_5;
	display y1_9_5;
	display z1_9_5;
	display x2_9_5;
	display y2_9_5;
	display z2_9_5;
	display x3_9_5;
	display y3_9_5;
	display z3_9_5;
	display x4_9_5;
	display y4_9_5;
	display z4_9_5;
	display x5_9_5;
	display y5_9_5;
	display z5_9_5;
	display x6_9_5;
	display y6_9_5;
	display z6_9_5;
	display x7_9_5;
	display y7_9_5;
	display z7_9_5;
	display x8_9_5;
	display y8_9_5;
	display z8_9_5;
	display x9_9_5;
	display y9_9_5;
	display z9_9_5;
	display x1_1_6;
	display y1_1_6;
	display z1_1_6;
	display x2_1_6;
	display y2_1_6;
	display z2_1_6;
	display x3_1_6;
	display y3_1_6;
	display z3_1_6;
	display x4_1_6;
	display y4_1_6;
	display z4_1_6;
	display x5_1_6;
	display y5_1_6;
	display z5_1_6;
	display x6_1_6;
	display y6_1_6;
	display z6_1_6;
	display x7_1_6;
	display y7_1_6;
	display z7_1_6;
	display x8_1_6;
	display y8_1_6;
	display z8_1_6;
	display x9_1_6;
	display y9_1_6;
	display z9_1_6;
	display x1_2_6;
	display y1_2_6;
	display z1_2_6;
	display x2_2_6;
	display y2_2_6;
	display z2_2_6;
	display x3_2_6;
	display y3_2_6;
	display z3_2_6;
	display x4_2_6;
	display y4_2_6;
	display z4_2_6;
	display x5_2_6;
	display y5_2_6;
	display z5_2_6;
	display x6_2_6;
	display y6_2_6;
	display z6_2_6;
	display x7_2_6;
	display y7_2_6;
	display z7_2_6;
	display x8_2_6;
	display y8_2_6;
	display z8_2_6;
	display x9_2_6;
	display y9_2_6;
	display z9_2_6;
	display x1_3_6;
	display y1_3_6;
	display z1_3_6;
	display x2_3_6;
	display y2_3_6;
	display z2_3_6;
	display x3_3_6;
	display y3_3_6;
	display z3_3_6;
	display x4_3_6;
	display y4_3_6;
	display z4_3_6;
	display x5_3_6;
	display y5_3_6;
	display z5_3_6;
	display x6_3_6;
	display y6_3_6;
	display z6_3_6;
	display x7_3_6;
	display y7_3_6;
	display z7_3_6;
	display x8_3_6;
	display y8_3_6;
	display z8_3_6;
	display x9_3_6;
	display y9_3_6;
	display z9_3_6;
	display x1_4_6;
	display y1_4_6;
	display z1_4_6;
	display x2_4_6;
	display y2_4_6;
	display z2_4_6;
	display x3_4_6;
	display y3_4_6;
	display z3_4_6;
	display x4_4_6;
	display y4_4_6;
	display z4_4_6;
	display x5_4_6;
	display y5_4_6;
	display z5_4_6;
	display x6_4_6;
	display y6_4_6;
	display z6_4_6;
	display x7_4_6;
	display y7_4_6;
	display z7_4_6;
	display x8_4_6;
	display y8_4_6;
	display z8_4_6;
	display x9_4_6;
	display y9_4_6;
	display z9_4_6;
	display x1_5_6;
	display y1_5_6;
	display z1_5_6;
	display x2_5_6;
	display y2_5_6;
	display z2_5_6;
	display x3_5_6;
	display y3_5_6;
	display z3_5_6;
	display x4_5_6;
	display y4_5_6;
	display z4_5_6;
	display x5_5_6;
	display y5_5_6;
	display z5_5_6;
	display x6_5_6;
	display y6_5_6;
	display z6_5_6;
	display x7_5_6;
	display y7_5_6;
	display z7_5_6;
	display x8_5_6;
	display y8_5_6;
	display z8_5_6;
	display x9_5_6;
	display y9_5_6;
	display z9_5_6;
	display x1_6_6;
	display y1_6_6;
	display z1_6_6;
	display x2_6_6;
	display y2_6_6;
	display z2_6_6;
	display x3_6_6;
	display y3_6_6;
	display z3_6_6;
	display x4_6_6;
	display y4_6_6;
	display z4_6_6;
	display x5_6_6;
	display y5_6_6;
	display z5_6_6;
	display x6_6_6;
	display y6_6_6;
	display z6_6_6;
	display x7_6_6;
	display y7_6_6;
	display z7_6_6;
	display x8_6_6;
	display y8_6_6;
	display z8_6_6;
	display x9_6_6;
	display y9_6_6;
	display z9_6_6;
	display x1_7_6;
	display y1_7_6;
	display z1_7_6;
	display x2_7_6;
	display y2_7_6;
	display z2_7_6;
	display x3_7_6;
	display y3_7_6;
	display z3_7_6;
	display x4_7_6;
	display y4_7_6;
	display z4_7_6;
	display x5_7_6;
	display y5_7_6;
	display z5_7_6;
	display x6_7_6;
	display y6_7_6;
	display z6_7_6;
	display x7_7_6;
	display y7_7_6;
	display z7_7_6;
	display x8_7_6;
	display y8_7_6;
	display z8_7_6;
	display x9_7_6;
	display y9_7_6;
	display z9_7_6;
	display x1_8_6;
	display y1_8_6;
	display z1_8_6;
	display x2_8_6;
	display y2_8_6;
	display z2_8_6;
	display x3_8_6;
	display y3_8_6;
	display z3_8_6;
	display x4_8_6;
	display y4_8_6;
	display z4_8_6;
	display x5_8_6;
	display y5_8_6;
	display z5_8_6;
	display x6_8_6;
	display y6_8_6;
	display z6_8_6;
	display x7_8_6;
	display y7_8_6;
	display z7_8_6;
	display x8_8_6;
	display y8_8_6;
	display z8_8_6;
	display x9_8_6;
	display y9_8_6;
	display z9_8_6;
	display x1_9_6;
	display y1_9_6;
	display z1_9_6;
	display x2_9_6;
	display y2_9_6;
	display z2_9_6;
	display x3_9_6;
	display y3_9_6;
	display z3_9_6;
	display x4_9_6;
	display y4_9_6;
	display z4_9_6;
	display x5_9_6;
	display y5_9_6;
	display z5_9_6;
	display x6_9_6;
	display y6_9_6;
	display z6_9_6;
	display x7_9_6;
	display y7_9_6;
	display z7_9_6;
	display x8_9_6;
	display y8_9_6;
	display z8_9_6;
	display x9_9_6;
	display y9_9_6;
	display z9_9_6;
	display x1_1_7;
	display y1_1_7;
	display z1_1_7;
	display x2_1_7;
	display y2_1_7;
	display z2_1_7;
	display x3_1_7;
	display y3_1_7;
	display z3_1_7;
	display x4_1_7;
	display y4_1_7;
	display z4_1_7;
	display x5_1_7;
	display y5_1_7;
	display z5_1_7;
	display x6_1_7;
	display y6_1_7;
	display z6_1_7;
	display x7_1_7;
	display y7_1_7;
	display z7_1_7;
	display x8_1_7;
	display y8_1_7;
	display z8_1_7;
	display x9_1_7;
	display y9_1_7;
	display z9_1_7;
	display x1_2_7;
	display y1_2_7;
	display z1_2_7;
	display x2_2_7;
	display y2_2_7;
	display z2_2_7;
	display x3_2_7;
	display y3_2_7;
	display z3_2_7;
	display x4_2_7;
	display y4_2_7;
	display z4_2_7;
	display x5_2_7;
	display y5_2_7;
	display z5_2_7;
	display x6_2_7;
	display y6_2_7;
	display z6_2_7;
	display x7_2_7;
	display y7_2_7;
	display z7_2_7;
	display x8_2_7;
	display y8_2_7;
	display z8_2_7;
	display x9_2_7;
	display y9_2_7;
	display z9_2_7;
	display x1_3_7;
	display y1_3_7;
	display z1_3_7;
	display x2_3_7;
	display y2_3_7;
	display z2_3_7;
	display x3_3_7;
	display y3_3_7;
	display z3_3_7;
	display x4_3_7;
	display y4_3_7;
	display z4_3_7;
	display x5_3_7;
	display y5_3_7;
	display z5_3_7;
	display x6_3_7;
	display y6_3_7;
	display z6_3_7;
	display x7_3_7;
	display y7_3_7;
	display z7_3_7;
	display x8_3_7;
	display y8_3_7;
	display z8_3_7;
	display x9_3_7;
	display y9_3_7;
	display z9_3_7;
	display x1_4_7;
	display y1_4_7;
	display z1_4_7;
	display x2_4_7;
	display y2_4_7;
	display z2_4_7;
	display x3_4_7;
	display y3_4_7;
	display z3_4_7;
	display x4_4_7;
	display y4_4_7;
	display z4_4_7;
	display x5_4_7;
	display y5_4_7;
	display z5_4_7;
	display x6_4_7;
	display y6_4_7;
	display z6_4_7;
	display x7_4_7;
	display y7_4_7;
	display z7_4_7;
	display x8_4_7;
	display y8_4_7;
	display z8_4_7;
	display x9_4_7;
	display y9_4_7;
	display z9_4_7;
	display x1_5_7;
	display y1_5_7;
	display z1_5_7;
	display x2_5_7;
	display y2_5_7;
	display z2_5_7;
	display x3_5_7;
	display y3_5_7;
	display z3_5_7;
	display x4_5_7;
	display y4_5_7;
	display z4_5_7;
	display x5_5_7;
	display y5_5_7;
	display z5_5_7;
	display x6_5_7;
	display y6_5_7;
	display z6_5_7;
	display x7_5_7;
	display y7_5_7;
	display z7_5_7;
	display x8_5_7;
	display y8_5_7;
	display z8_5_7;
	display x9_5_7;
	display y9_5_7;
	display z9_5_7;
	display x1_6_7;
	display y1_6_7;
	display z1_6_7;
	display x2_6_7;
	display y2_6_7;
	display z2_6_7;
	display x3_6_7;
	display y3_6_7;
	display z3_6_7;
	display x4_6_7;
	display y4_6_7;
	display z4_6_7;
	display x5_6_7;
	display y5_6_7;
	display z5_6_7;
	display x6_6_7;
	display y6_6_7;
	display z6_6_7;
	display x7_6_7;
	display y7_6_7;
	display z7_6_7;
	display x8_6_7;
	display y8_6_7;
	display z8_6_7;
	display x9_6_7;
	display y9_6_7;
	display z9_6_7;
	display x1_7_7;
	display y1_7_7;
	display z1_7_7;
	display x2_7_7;
	display y2_7_7;
	display z2_7_7;
	display x3_7_7;
	display y3_7_7;
	display z3_7_7;
	display x4_7_7;
	display y4_7_7;
	display z4_7_7;
	display x5_7_7;
	display y5_7_7;
	display z5_7_7;
	display x6_7_7;
	display y6_7_7;
	display z6_7_7;
	display x7_7_7;
	display y7_7_7;
	display z7_7_7;
	display x8_7_7;
	display y8_7_7;
	display z8_7_7;
	display x9_7_7;
	display y9_7_7;
	display z9_7_7;
	display x1_8_7;
	display y1_8_7;
	display z1_8_7;
	display x2_8_7;
	display y2_8_7;
	display z2_8_7;
	display x3_8_7;
	display y3_8_7;
	display z3_8_7;
	display x4_8_7;
	display y4_8_7;
	display z4_8_7;
	display x5_8_7;
	display y5_8_7;
	display z5_8_7;
	display x6_8_7;
	display y6_8_7;
	display z6_8_7;
	display x7_8_7;
	display y7_8_7;
	display z7_8_7;
	display x8_8_7;
	display y8_8_7;
	display z8_8_7;
	display x9_8_7;
	display y9_8_7;
	display z9_8_7;
	display x1_9_7;
	display y1_9_7;
	display z1_9_7;
	display x2_9_7;
	display y2_9_7;
	display z2_9_7;
	display x3_9_7;
	display y3_9_7;
	display z3_9_7;
	display x4_9_7;
	display y4_9_7;
	display z4_9_7;
	display x5_9_7;
	display y5_9_7;
	display z5_9_7;
	display x6_9_7;
	display y6_9_7;
	display z6_9_7;
	display x7_9_7;
	display y7_9_7;
	display z7_9_7;
	display x8_9_7;
	display y8_9_7;
	display z8_9_7;
	display x9_9_7;
	display y9_9_7;
	display z9_9_7;
	display x1_1_8;
	display y1_1_8;
	display z1_1_8;
	display x2_1_8;
	display y2_1_8;
	display z2_1_8;
	display x3_1_8;
	display y3_1_8;
	display z3_1_8;
	display x4_1_8;
	display y4_1_8;
	display z4_1_8;
	display x5_1_8;
	display y5_1_8;
	display z5_1_8;
	display x6_1_8;
	display y6_1_8;
	display z6_1_8;
	display x7_1_8;
	display y7_1_8;
	display z7_1_8;
	display x8_1_8;
	display y8_1_8;
	display z8_1_8;
	display x9_1_8;
	display y9_1_8;
	display z9_1_8;
	display x1_2_8;
	display y1_2_8;
	display z1_2_8;
	display x2_2_8;
	display y2_2_8;
	display z2_2_8;
	display x3_2_8;
	display y3_2_8;
	display z3_2_8;
	display x4_2_8;
	display y4_2_8;
	display z4_2_8;
	display x5_2_8;
	display y5_2_8;
	display z5_2_8;
	display x6_2_8;
	display y6_2_8;
	display z6_2_8;
	display x7_2_8;
	display y7_2_8;
	display z7_2_8;
	display x8_2_8;
	display y8_2_8;
	display z8_2_8;
	display x9_2_8;
	display y9_2_8;
	display z9_2_8;
	display x1_3_8;
	display y1_3_8;
	display z1_3_8;
	display x2_3_8;
	display y2_3_8;
	display z2_3_8;
	display x3_3_8;
	display y3_3_8;
	display z3_3_8;
	display x4_3_8;
	display y4_3_8;
	display z4_3_8;
	display x5_3_8;
	display y5_3_8;
	display z5_3_8;
	display x6_3_8;
	display y6_3_8;
	display z6_3_8;
	display x7_3_8;
	display y7_3_8;
	display z7_3_8;
	display x8_3_8;
	display y8_3_8;
	display z8_3_8;
	display x9_3_8;
	display y9_3_8;
	display z9_3_8;
	display x1_4_8;
	display y1_4_8;
	display z1_4_8;
	display x2_4_8;
	display y2_4_8;
	display z2_4_8;
	display x3_4_8;
	display y3_4_8;
	display z3_4_8;
	display x4_4_8;
	display y4_4_8;
	display z4_4_8;
	display x5_4_8;
	display y5_4_8;
	display z5_4_8;
	display x6_4_8;
	display y6_4_8;
	display z6_4_8;
	display x7_4_8;
	display y7_4_8;
	display z7_4_8;
	display x8_4_8;
	display y8_4_8;
	display z8_4_8;
	display x9_4_8;
	display y9_4_8;
	display z9_4_8;
	display x1_5_8;
	display y1_5_8;
	display z1_5_8;
	display x2_5_8;
	display y2_5_8;
	display z2_5_8;
	display x3_5_8;
	display y3_5_8;
	display z3_5_8;
	display x4_5_8;
	display y4_5_8;
	display z4_5_8;
	display x5_5_8;
	display y5_5_8;
	display z5_5_8;
	display x6_5_8;
	display y6_5_8;
	display z6_5_8;
	display x7_5_8;
	display y7_5_8;
	display z7_5_8;
	display x8_5_8;
	display y8_5_8;
	display z8_5_8;
	display x9_5_8;
	display y9_5_8;
	display z9_5_8;
	display x1_6_8;
	display y1_6_8;
	display z1_6_8;
	display x2_6_8;
	display y2_6_8;
	display z2_6_8;
	display x3_6_8;
	display y3_6_8;
	display z3_6_8;
	display x4_6_8;
	display y4_6_8;
	display z4_6_8;
	display x5_6_8;
	display y5_6_8;
	display z5_6_8;
	display x6_6_8;
	display y6_6_8;
	display z6_6_8;
	display x7_6_8;
	display y7_6_8;
	display z7_6_8;
	display x8_6_8;
	display y8_6_8;
	display z8_6_8;
	display x9_6_8;
	display y9_6_8;
	display z9_6_8;
	display x1_7_8;
	display y1_7_8;
	display z1_7_8;
	display x2_7_8;
	display y2_7_8;
	display z2_7_8;
	display x3_7_8;
	display y3_7_8;
	display z3_7_8;
	display x4_7_8;
	display y4_7_8;
	display z4_7_8;
	display x5_7_8;
	display y5_7_8;
	display z5_7_8;
	display x6_7_8;
	display y6_7_8;
	display z6_7_8;
	display x7_7_8;
	display y7_7_8;
	display z7_7_8;
	display x8_7_8;
	display y8_7_8;
	display z8_7_8;
	display x9_7_8;
	display y9_7_8;
	display z9_7_8;
	display x1_8_8;
	display y1_8_8;
	display z1_8_8;
	display x2_8_8;
	display y2_8_8;
	display z2_8_8;
	display x3_8_8;
	display y3_8_8;
	display z3_8_8;
	display x4_8_8;
	display y4_8_8;
	display z4_8_8;
	display x5_8_8;
	display y5_8_8;
	display z5_8_8;
	display x6_8_8;
	display y6_8_8;
	display z6_8_8;
	display x7_8_8;
	display y7_8_8;
	display z7_8_8;
	display x8_8_8;
	display y8_8_8;
	display z8_8_8;
	display x9_8_8;
	display y9_8_8;
	display z9_8_8;
	display x1_9_8;
	display y1_9_8;
	display z1_9_8;
	display x2_9_8;
	display y2_9_8;
	display z2_9_8;
	display x3_9_8;
	display y3_9_8;
	display z3_9_8;
	display x4_9_8;
	display y4_9_8;
	display z4_9_8;
	display x5_9_8;
	display y5_9_8;
	display z5_9_8;
	display x6_9_8;
	display y6_9_8;
	display z6_9_8;
	display x7_9_8;
	display y7_9_8;
	display z7_9_8;
	display x8_9_8;
	display y8_9_8;
	display z8_9_8;
	display x9_9_8;
	display y9_9_8;
	display z9_9_8;
	display x1_1_9;
	display y1_1_9;
	display z1_1_9;
	display x2_1_9;
	display y2_1_9;
	display z2_1_9;
	display x3_1_9;
	display y3_1_9;
	display z3_1_9;
	display x4_1_9;
	display y4_1_9;
	display z4_1_9;
	display x5_1_9;
	display y5_1_9;
	display z5_1_9;
	display x6_1_9;
	display y6_1_9;
	display z6_1_9;
	display x7_1_9;
	display y7_1_9;
	display z7_1_9;
	display x8_1_9;
	display y8_1_9;
	display z8_1_9;
	display x9_1_9;
	display y9_1_9;
	display z9_1_9;
	display x1_2_9;
	display y1_2_9;
	display z1_2_9;
	display x2_2_9;
	display y2_2_9;
	display z2_2_9;
	display x3_2_9;
	display y3_2_9;
	display z3_2_9;
	display x4_2_9;
	display y4_2_9;
	display z4_2_9;
	display x5_2_9;
	display y5_2_9;
	display z5_2_9;
	display x6_2_9;
	display y6_2_9;
	display z6_2_9;
	display x7_2_9;
	display y7_2_9;
	display z7_2_9;
	display x8_2_9;
	display y8_2_9;
	display z8_2_9;
	display x9_2_9;
	display y9_2_9;
	display z9_2_9;
	display x1_3_9;
	display y1_3_9;
	display z1_3_9;
	display x2_3_9;
	display y2_3_9;
	display z2_3_9;
	display x3_3_9;
	display y3_3_9;
	display z3_3_9;
	display x4_3_9;
	display y4_3_9;
	display z4_3_9;
	display x5_3_9;
	display y5_3_9;
	display z5_3_9;
	display x6_3_9;
	display y6_3_9;
	display z6_3_9;
	display x7_3_9;
	display y7_3_9;
	display z7_3_9;
	display x8_3_9;
	display y8_3_9;
	display z8_3_9;
	display x9_3_9;
	display y9_3_9;
	display z9_3_9;
	display x1_4_9;
	display y1_4_9;
	display z1_4_9;
	display x2_4_9;
	display y2_4_9;
	display z2_4_9;
	display x3_4_9;
	display y3_4_9;
	display z3_4_9;
	display x4_4_9;
	display y4_4_9;
	display z4_4_9;
	display x5_4_9;
	display y5_4_9;
	display z5_4_9;
	display x6_4_9;
	display y6_4_9;
	display z6_4_9;
	display x7_4_9;
	display y7_4_9;
	display z7_4_9;
	display x8_4_9;
	display y8_4_9;
	display z8_4_9;
	display x9_4_9;
	display y9_4_9;
	display z9_4_9;
	display x1_5_9;
	display y1_5_9;
	display z1_5_9;
	display x2_5_9;
	display y2_5_9;
	display z2_5_9;
	display x3_5_9;
	display y3_5_9;
	display z3_5_9;
	display x4_5_9;
	display y4_5_9;
	display z4_5_9;
	display x5_5_9;
	display y5_5_9;
	display z5_5_9;
	display x6_5_9;
	display y6_5_9;
	display z6_5_9;
	display x7_5_9;
	display y7_5_9;
	display z7_5_9;
	display x8_5_9;
	display y8_5_9;
	display z8_5_9;
	display x9_5_9;
	display y9_5_9;
	display z9_5_9;
	display x1_6_9;
	display y1_6_9;
	display z1_6_9;
	display x2_6_9;
	display y2_6_9;
	display z2_6_9;
	display x3_6_9;
	display y3_6_9;
	display z3_6_9;
	display x4_6_9;
	display y4_6_9;
	display z4_6_9;
	display x5_6_9;
	display y5_6_9;
	display z5_6_9;
	display x6_6_9;
	display y6_6_9;
	display z6_6_9;
	display x7_6_9;
	display y7_6_9;
	display z7_6_9;
	display x8_6_9;
	display y8_6_9;
	display z8_6_9;
	display x9_6_9;
	display y9_6_9;
	display z9_6_9;
	display x1_7_9;
	display y1_7_9;
	display z1_7_9;
	display x2_7_9;
	display y2_7_9;
	display z2_7_9;
	display x3_7_9;
	display y3_7_9;
	display z3_7_9;
	display x4_7_9;
	display y4_7_9;
	display z4_7_9;
	display x5_7_9;
	display y5_7_9;
	display z5_7_9;
	display x6_7_9;
	display y6_7_9;
	display z6_7_9;
	display x7_7_9;
	display y7_7_9;
	display z7_7_9;
	display x8_7_9;
	display y8_7_9;
	display z8_7_9;
	display x9_7_9;
	display y9_7_9;
	display z9_7_9;
	display x1_8_9;
	display y1_8_9;
	display z1_8_9;
	display x2_8_9;
	display y2_8_9;
	display z2_8_9;
	display x3_8_9;
	display y3_8_9;
	display z3_8_9;
	display x4_8_9;
	display y4_8_9;
	display z4_8_9;
	display x5_8_9;
	display y5_8_9;
	display z5_8_9;
	display x6_8_9;
	display y6_8_9;
	display z6_8_9;
	display x7_8_9;
	display y7_8_9;
	display z7_8_9;
	display x8_8_9;
	display y8_8_9;
	display z8_8_9;
	display x9_8_9;
	display y9_8_9;
	display z9_8_9;
	display x1_9_9;
	display y1_9_9;
	display z1_9_9;
	display x2_9_9;
	display y2_9_9;
	display z2_9_9;
	display x3_9_9;
	display y3_9_9;
	display z3_9_9;
	display x4_9_9;
	display y4_9_9;
	display z4_9_9;
	display x5_9_9;
	display y5_9_9;
	display z5_9_9;
	display x6_9_9;
	display y6_9_9;
	display z6_9_9;
	display x7_9_9;
	display y7_9_9;
	display z7_9_9;
	display x8_9_9;
	display y8_9_9;
	display z8_9_9;
	display x9_9_9;
	display y9_9_9;
	display z9_9_9;
	display y10_1_1;
	display z10_1_1;
	display y10_2_1;
	display z10_2_1;
	display y10_3_1;
	display z10_3_1;
	display y10_4_1;
	display z10_4_1;
	display y10_5_1;
	display z10_5_1;
	display y10_6_1;
	display z10_6_1;
	display y10_7_1;
	display z10_7_1;
	display y10_8_1;
	display z10_8_1;
	display y10_9_1;
	display z10_9_1;
	display y10_1_2;
	display z10_1_2;
	display y10_2_2;
	display z10_2_2;
	display y10_3_2;
	display z10_3_2;
	display y10_4_2;
	display z10_4_2;
	display y10_5_2;
	display z10_5_2;
	display y10_6_2;
	display z10_6_2;
	display y10_7_2;
	display z10_7_2;
	display y10_8_2;
	display z10_8_2;
	display y10_9_2;
	display z10_9_2;
	display y10_1_3;
	display z10_1_3;
	display y10_2_3;
	display z10_2_3;
	display y10_3_3;
	display z10_3_3;
	display y10_4_3;
	display z10_4_3;
	display y10_5_3;
	display z10_5_3;
	display y10_6_3;
	display z10_6_3;
	display y10_7_3;
	display z10_7_3;
	display y10_8_3;
	display z10_8_3;
	display y10_9_3;
	display z10_9_3;
	display y10_1_4;
	display z10_1_4;
	display y10_2_4;
	display z10_2_4;
	display y10_3_4;
	display z10_3_4;
	display y10_4_4;
	display z10_4_4;
	display y10_5_4;
	display z10_5_4;
	display y10_6_4;
	display z10_6_4;
	display y10_7_4;
	display z10_7_4;
	display y10_8_4;
	display z10_8_4;
	display y10_9_4;
	display z10_9_4;
	display y10_1_5;
	display z10_1_5;
	display y10_2_5;
	display z10_2_5;
	display y10_3_5;
	display z10_3_5;
	display y10_4_5;
	display z10_4_5;
	display y10_5_5;
	display z10_5_5;
	display y10_6_5;
	display z10_6_5;
	display y10_7_5;
	display z10_7_5;
	display y10_8_5;
	display z10_8_5;
	display y10_9_5;
	display z10_9_5;
	display y10_1_6;
	display z10_1_6;
	display y10_2_6;
	display z10_2_6;
	display y10_3_6;
	display z10_3_6;
	display y10_4_6;
	display z10_4_6;
	display y10_5_6;
	display z10_5_6;
	display y10_6_6;
	display z10_6_6;
	display y10_7_6;
	display z10_7_6;
	display y10_8_6;
	display z10_8_6;
	display y10_9_6;
	display z10_9_6;
	display y10_1_7;
	display z10_1_7;
	display y10_2_7;
	display z10_2_7;
	display y10_3_7;
	display z10_3_7;
	display y10_4_7;
	display z10_4_7;
	display y10_5_7;
	display z10_5_7;
	display y10_6_7;
	display z10_6_7;
	display y10_7_7;
	display z10_7_7;
	display y10_8_7;
	display z10_8_7;
	display y10_9_7;
	display z10_9_7;
	display y10_1_8;
	display z10_1_8;
	display y10_2_8;
	display z10_2_8;
	display y10_3_8;
	display z10_3_8;
	display y10_4_8;
	display z10_4_8;
	display y10_5_8;
	display z10_5_8;
	display y10_6_8;
	display z10_6_8;
	display y10_7_8;
	display z10_7_8;
	display y10_8_8;
	display z10_8_8;
	display y10_9_8;
	display z10_9_8;
	display y10_1_9;
	display z10_1_9;
	display y10_2_9;
	display z10_2_9;
	display y10_3_9;
	display z10_3_9;
	display y10_4_9;
	display z10_4_9;
	display y10_5_9;
	display z10_5_9;
	display y10_6_9;
	display z10_6_9;
	display y10_7_9;
	display z10_7_9;
	display y10_8_9;
	display z10_8_9;
	display y10_9_9;
	display z10_9_9;
	display x1_10_1;
	display z1_10_1;
	display x2_10_1;
	display z2_10_1;
	display x3_10_1;
	display z3_10_1;
	display x4_10_1;
	display z4_10_1;
	display x5_10_1;
	display z5_10_1;
	display x6_10_1;
	display z6_10_1;
	display x7_10_1;
	display z7_10_1;
	display x8_10_1;
	display z8_10_1;
	display x9_10_1;
	display z9_10_1;
	display x1_10_2;
	display z1_10_2;
	display x2_10_2;
	display z2_10_2;
	display x3_10_2;
	display z3_10_2;
	display x4_10_2;
	display z4_10_2;
	display x5_10_2;
	display z5_10_2;
	display x6_10_2;
	display z6_10_2;
	display x7_10_2;
	display z7_10_2;
	display x8_10_2;
	display z8_10_2;
	display x9_10_2;
	display z9_10_2;
	display x1_10_3;
	display z1_10_3;
	display x2_10_3;
	display z2_10_3;
	display x3_10_3;
	display z3_10_3;
	display x4_10_3;
	display z4_10_3;
	display x5_10_3;
	display z5_10_3;
	display x6_10_3;
	display z6_10_3;
	display x7_10_3;
	display z7_10_3;
	display x8_10_3;
	display z8_10_3;
	display x9_10_3;
	display z9_10_3;
	display x1_10_4;
	display z1_10_4;
	display x2_10_4;
	display z2_10_4;
	display x3_10_4;
	display z3_10_4;
	display x4_10_4;
	display z4_10_4;
	display x5_10_4;
	display z5_10_4;
	display x6_10_4;
	display z6_10_4;
	display x7_10_4;
	display z7_10_4;
	display x8_10_4;
	display z8_10_4;
	display x9_10_4;
	display z9_10_4;
	display x1_10_5;
	display z1_10_5;
	display x2_10_5;
	display z2_10_5;
	display x3_10_5;
	display z3_10_5;
	display x4_10_5;
	display z4_10_5;
	display x5_10_5;
	display z5_10_5;
	display x6_10_5;
	display z6_10_5;
	display x7_10_5;
	display z7_10_5;
	display x8_10_5;
	display z8_10_5;
	display x9_10_5;
	display z9_10_5;
	display x1_10_6;
	display z1_10_6;
	display x2_10_6;
	display z2_10_6;
	display x3_10_6;
	display z3_10_6;
	display x4_10_6;
	display z4_10_6;
	display x5_10_6;
	display z5_10_6;
	display x6_10_6;
	display z6_10_6;
	display x7_10_6;
	display z7_10_6;
	display x8_10_6;
	display z8_10_6;
	display x9_10_6;
	display z9_10_6;
	display x1_10_7;
	display z1_10_7;
	display x2_10_7;
	display z2_10_7;
	display x3_10_7;
	display z3_10_7;
	display x4_10_7;
	display z4_10_7;
	display x5_10_7;
	display z5_10_7;
	display x6_10_7;
	display z6_10_7;
	display x7_10_7;
	display z7_10_7;
	display x8_10_7;
	display z8_10_7;
	display x9_10_7;
	display z9_10_7;
	display x1_10_8;
	display z1_10_8;
	display x2_10_8;
	display z2_10_8;
	display x3_10_8;
	display z3_10_8;
	display x4_10_8;
	display z4_10_8;
	display x5_10_8;
	display z5_10_8;
	display x6_10_8;
	display z6_10_8;
	display x7_10_8;
	display z7_10_8;
	display x8_10_8;
	display z8_10_8;
	display x9_10_8;
	display z9_10_8;
	display x1_10_9;
	display z1_10_9;
	display x2_10_9;
	display z2_10_9;
	display x3_10_9;
	display z3_10_9;
	display x4_10_9;
	display z4_10_9;
	display x5_10_9;
	display z5_10_9;
	display x6_10_9;
	display z6_10_9;
	display x7_10_9;
	display z7_10_9;
	display x8_10_9;
	display z8_10_9;
	display x9_10_9;
	display z9_10_9;
	display x1_1_10;
	display y1_1_10;
	display x2_1_10;
	display y2_1_10;
	display x3_1_10;
	display y3_1_10;
	display x4_1_10;
	display y4_1_10;
	display x5_1_10;
	display y5_1_10;
	display x6_1_10;
	display y6_1_10;
	display x7_1_10;
	display y7_1_10;
	display x8_1_10;
	display y8_1_10;
	display x9_1_10;
	display y9_1_10;
	display x1_2_10;
	display y1_2_10;
	display x2_2_10;
	display y2_2_10;
	display x3_2_10;
	display y3_2_10;
	display x4_2_10;
	display y4_2_10;
	display x5_2_10;
	display y5_2_10;
	display x6_2_10;
	display y6_2_10;
	display x7_2_10;
	display y7_2_10;
	display x8_2_10;
	display y8_2_10;
	display x9_2_10;
	display y9_2_10;
	display x1_3_10;
	display y1_3_10;
	display x2_3_10;
	display y2_3_10;
	display x3_3_10;
	display y3_3_10;
	display x4_3_10;
	display y4_3_10;
	display x5_3_10;
	display y5_3_10;
	display x6_3_10;
	display y6_3_10;
	display x7_3_10;
	display y7_3_10;
	display x8_3_10;
	display y8_3_10;
	display x9_3_10;
	display y9_3_10;
	display x1_4_10;
	display y1_4_10;
	display x2_4_10;
	display y2_4_10;
	display x3_4_10;
	display y3_4_10;
	display x4_4_10;
	display y4_4_10;
	display x5_4_10;
	display y5_4_10;
	display x6_4_10;
	display y6_4_10;
	display x7_4_10;
	display y7_4_10;
	display x8_4_10;
	display y8_4_10;
	display x9_4_10;
	display y9_4_10;
	display x1_5_10;
	display y1_5_10;
	display x2_5_10;
	display y2_5_10;
	display x3_5_10;
	display y3_5_10;
	display x4_5_10;
	display y4_5_10;
	display x5_5_10;
	display y5_5_10;
	display x6_5_10;
	display y6_5_10;
	display x7_5_10;
	display y7_5_10;
	display x8_5_10;
	display y8_5_10;
	display x9_5_10;
	display y9_5_10;
	display x1_6_10;
	display y1_6_10;
	display x2_6_10;
	display y2_6_10;
	display x3_6_10;
	display y3_6_10;
	display x4_6_10;
	display y4_6_10;
	display x5_6_10;
	display y5_6_10;
	display x6_6_10;
	display y6_6_10;
	display x7_6_10;
	display y7_6_10;
	display x8_6_10;
	display y8_6_10;
	display x9_6_10;
	display y9_6_10;
	display x1_7_10;
	display y1_7_10;
	display x2_7_10;
	display y2_7_10;
	display x3_7_10;
	display y3_7_10;
	display x4_7_10;
	display y4_7_10;
	display x5_7_10;
	display y5_7_10;
	display x6_7_10;
	display y6_7_10;
	display x7_7_10;
	display y7_7_10;
	display x8_7_10;
	display y8_7_10;
	display x9_7_10;
	display y9_7_10;
	display x1_8_10;
	display y1_8_10;
	display x2_8_10;
	display y2_8_10;
	display x3_8_10;
	display y3_8_10;
	display x4_8_10;
	display y4_8_10;
	display x5_8_10;
	display y5_8_10;
	display x6_8_10;
	display y6_8_10;
	display x7_8_10;
	display y7_8_10;
	display x8_8_10;
	display y8_8_10;
	display x9_8_10;
	display y9_8_10;
	display x1_9_10;
	display y1_9_10;
	display x2_9_10;
	display y2_9_10;
	display x3_9_10;
	display y3_9_10;
	display x4_9_10;
	display y4_9_10;
	display x5_9_10;
	display y5_9_10;
	display x6_9_10;
	display y6_9_10;
	display x7_9_10;
	display y7_9_10;
	display x8_9_10;
	display y8_9_10;
	display x9_9_10;
	display y9_9_10;
	display y0_1_1;
	display z0_1_1;
	display y11_1_1;
	display z11_1_1;
	display y0_2_1;
	display z0_2_1;
	display y11_2_1;
	display z11_2_1;
	display y0_3_1;
	display z0_3_1;
	display y11_3_1;
	display z11_3_1;
	display y0_4_1;
	display z0_4_1;
	display y11_4_1;
	display z11_4_1;
	display y0_5_1;
	display z0_5_1;
	display y11_5_1;
	display z11_5_1;
	display y0_6_1;
	display z0_6_1;
	display y11_6_1;
	display z11_6_1;
	display y0_7_1;
	display z0_7_1;
	display y11_7_1;
	display z11_7_1;
	display y0_8_1;
	display z0_8_1;
	display y11_8_1;
	display z11_8_1;
	display y0_9_1;
	display z0_9_1;
	display y11_9_1;
	display z11_9_1;
	display y0_10_1;
	display z0_10_1;
	display y11_10_1;
	display z11_10_1;
	display y0_1_2;
	display z0_1_2;
	display y11_1_2;
	display z11_1_2;
	display y0_2_2;
	display z0_2_2;
	display y11_2_2;
	display z11_2_2;
	display y0_3_2;
	display z0_3_2;
	display y11_3_2;
	display z11_3_2;
	display y0_4_2;
	display z0_4_2;
	display y11_4_2;
	display z11_4_2;
	display y0_5_2;
	display z0_5_2;
	display y11_5_2;
	display z11_5_2;
	display y0_6_2;
	display z0_6_2;
	display y11_6_2;
	display z11_6_2;
	display y0_7_2;
	display z0_7_2;
	display y11_7_2;
	display z11_7_2;
	display y0_8_2;
	display z0_8_2;
	display y11_8_2;
	display z11_8_2;
	display y0_9_2;
	display z0_9_2;
	display y11_9_2;
	display z11_9_2;
	display y0_10_2;
	display z0_10_2;
	display y11_10_2;
	display z11_10_2;
	display y0_1_3;
	display z0_1_3;
	display y11_1_3;
	display z11_1_3;
	display y0_2_3;
	display z0_2_3;
	display y11_2_3;
	display z11_2_3;
	display y0_3_3;
	display z0_3_3;
	display y11_3_3;
	display z11_3_3;
	display y0_4_3;
	display z0_4_3;
	display y11_4_3;
	display z11_4_3;
	display y0_5_3;
	display z0_5_3;
	display y11_5_3;
	display z11_5_3;
	display y0_6_3;
	display z0_6_3;
	display y11_6_3;
	display z11_6_3;
	display y0_7_3;
	display z0_7_3;
	display y11_7_3;
	display z11_7_3;
	display y0_8_3;
	display z0_8_3;
	display y11_8_3;
	display z11_8_3;
	display y0_9_3;
	display z0_9_3;
	display y11_9_3;
	display z11_9_3;
	display y0_10_3;
	display z0_10_3;
	display y11_10_3;
	display z11_10_3;
	display y0_1_4;
	display z0_1_4;
	display y11_1_4;
	display z11_1_4;
	display y0_2_4;
	display z0_2_4;
	display y11_2_4;
	display z11_2_4;
	display y0_3_4;
	display z0_3_4;
	display y11_3_4;
	display z11_3_4;
	display y0_4_4;
	display z0_4_4;
	display y11_4_4;
	display z11_4_4;
	display y0_5_4;
	display z0_5_4;
	display y11_5_4;
	display z11_5_4;
	display y0_6_4;
	display z0_6_4;
	display y11_6_4;
	display z11_6_4;
	display y0_7_4;
	display z0_7_4;
	display y11_7_4;
	display z11_7_4;
	display y0_8_4;
	display z0_8_4;
	display y11_8_4;
	display z11_8_4;
	display y0_9_4;
	display z0_9_4;
	display y11_9_4;
	display z11_9_4;
	display y0_10_4;
	display z0_10_4;
	display y11_10_4;
	display z11_10_4;
	display y0_1_5;
	display z0_1_5;
	display y11_1_5;
	display z11_1_5;
	display y0_2_5;
	display z0_2_5;
	display y11_2_5;
	display z11_2_5;
	display y0_3_5;
	display z0_3_5;
	display y11_3_5;
	display z11_3_5;
	display y0_4_5;
	display z0_4_5;
	display y11_4_5;
	display z11_4_5;
	display y0_5_5;
	display z0_5_5;
	display y11_5_5;
	display z11_5_5;
	display y0_6_5;
	display z0_6_5;
	display y11_6_5;
	display z11_6_5;
	display y0_7_5;
	display z0_7_5;
	display y11_7_5;
	display z11_7_5;
	display y0_8_5;
	display z0_8_5;
	display y11_8_5;
	display z11_8_5;
	display y0_9_5;
	display z0_9_5;
	display y11_9_5;
	display z11_9_5;
	display y0_10_5;
	display z0_10_5;
	display y11_10_5;
	display z11_10_5;
	display y0_1_6;
	display z0_1_6;
	display y11_1_6;
	display z11_1_6;
	display y0_2_6;
	display z0_2_6;
	display y11_2_6;
	display z11_2_6;
	display y0_3_6;
	display z0_3_6;
	display y11_3_6;
	display z11_3_6;
	display y0_4_6;
	display z0_4_6;
	display y11_4_6;
	display z11_4_6;
	display y0_5_6;
	display z0_5_6;
	display y11_5_6;
	display z11_5_6;
	display y0_6_6;
	display z0_6_6;
	display y11_6_6;
	display z11_6_6;
	display y0_7_6;
	display z0_7_6;
	display y11_7_6;
	display z11_7_6;
	display y0_8_6;
	display z0_8_6;
	display y11_8_6;
	display z11_8_6;
	display y0_9_6;
	display z0_9_6;
	display y11_9_6;
	display z11_9_6;
	display y0_10_6;
	display z0_10_6;
	display y11_10_6;
	display z11_10_6;
	display y0_1_7;
	display z0_1_7;
	display y11_1_7;
	display z11_1_7;
	display y0_2_7;
	display z0_2_7;
	display y11_2_7;
	display z11_2_7;
	display y0_3_7;
	display z0_3_7;
	display y11_3_7;
	display z11_3_7;
	display y0_4_7;
	display z0_4_7;
	display y11_4_7;
	display z11_4_7;
	display y0_5_7;
	display z0_5_7;
	display y11_5_7;
	display z11_5_7;
	display y0_6_7;
	display z0_6_7;
	display y11_6_7;
	display z11_6_7;
	display y0_7_7;
	display z0_7_7;
	display y11_7_7;
	display z11_7_7;
	display y0_8_7;
	display z0_8_7;
	display y11_8_7;
	display z11_8_7;
	display y0_9_7;
	display z0_9_7;
	display y11_9_7;
	display z11_9_7;
	display y0_10_7;
	display z0_10_7;
	display y11_10_7;
	display z11_10_7;
	display y0_1_8;
	display z0_1_8;
	display y11_1_8;
	display z11_1_8;
	display y0_2_8;
	display z0_2_8;
	display y11_2_8;
	display z11_2_8;
	display y0_3_8;
	display z0_3_8;
	display y11_3_8;
	display z11_3_8;
	display y0_4_8;
	display z0_4_8;
	display y11_4_8;
	display z11_4_8;
	display y0_5_8;
	display z0_5_8;
	display y11_5_8;
	display z11_5_8;
	display y0_6_8;
	display z0_6_8;
	display y11_6_8;
	display z11_6_8;
	display y0_7_8;
	display z0_7_8;
	display y11_7_8;
	display z11_7_8;
	display y0_8_8;
	display z0_8_8;
	display y11_8_8;
	display z11_8_8;
	display y0_9_8;
	display z0_9_8;
	display y11_9_8;
	display z11_9_8;
	display y0_10_8;
	display z0_10_8;
	display y11_10_8;
	display z11_10_8;
	display y0_1_9;
	display z0_1_9;
	display y11_1_9;
	display z11_1_9;
	display y0_2_9;
	display z0_2_9;
	display y11_2_9;
	display z11_2_9;
	display y0_3_9;
	display z0_3_9;
	display y11_3_9;
	display z11_3_9;
	display y0_4_9;
	display z0_4_9;
	display y11_4_9;
	display z11_4_9;
	display y0_5_9;
	display z0_5_9;
	display y11_5_9;
	display z11_5_9;
	display y0_6_9;
	display z0_6_9;
	display y11_6_9;
	display z11_6_9;
	display y0_7_9;
	display z0_7_9;
	display y11_7_9;
	display z11_7_9;
	display y0_8_9;
	display z0_8_9;
	display y11_8_9;
	display z11_8_9;
	display y0_9_9;
	display z0_9_9;
	display y11_9_9;
	display z11_9_9;
	display y0_10_9;
	display z0_10_9;
	display y11_10_9;
	display z11_10_9;
	display y0_1_10;
	display z0_1_10;
	display y11_1_10;
	display z11_1_10;
	display y0_2_10;
	display z0_2_10;
	display y11_2_10;
	display z11_2_10;
	display y0_3_10;
	display z0_3_10;
	display y11_3_10;
	display z11_3_10;
	display y0_4_10;
	display z0_4_10;
	display y11_4_10;
	display z11_4_10;
	display y0_5_10;
	display z0_5_10;
	display y11_5_10;
	display z11_5_10;
	display y0_6_10;
	display z0_6_10;
	display y11_6_10;
	display z11_6_10;
	display y0_7_10;
	display z0_7_10;
	display y11_7_10;
	display z11_7_10;
	display y0_8_10;
	display z0_8_10;
	display y11_8_10;
	display z11_8_10;
	display y0_9_10;
	display z0_9_10;
	display y11_9_10;
	display z11_9_10;
	display y0_10_10;
	display z0_10_10;
	display y11_10_10;
	display z11_10_10;
	display x1_0_1;
	display z1_0_1;
	display x1_11_1;
	display z1_11_1;
	display x2_0_1;
	display z2_0_1;
	display x2_11_1;
	display z2_11_1;
	display x3_0_1;
	display z3_0_1;
	display x3_11_1;
	display z3_11_1;
	display x4_0_1;
	display z4_0_1;
	display x4_11_1;
	display z4_11_1;
	display x5_0_1;
	display z5_0_1;
	display x5_11_1;
	display z5_11_1;
	display x6_0_1;
	display z6_0_1;
	display x6_11_1;
	display z6_11_1;
	display x7_0_1;
	display z7_0_1;
	display x7_11_1;
	display z7_11_1;
	display x8_0_1;
	display z8_0_1;
	display x8_11_1;
	display z8_11_1;
	display x9_0_1;
	display z9_0_1;
	display x9_11_1;
	display z9_11_1;
	display x10_0_1;
	display z10_0_1;
	display x10_11_1;
	display z10_11_1;
	display x1_0_2;
	display z1_0_2;
	display x1_11_2;
	display z1_11_2;
	display x2_0_2;
	display z2_0_2;
	display x2_11_2;
	display z2_11_2;
	display x3_0_2;
	display z3_0_2;
	display x3_11_2;
	display z3_11_2;
	display x4_0_2;
	display z4_0_2;
	display x4_11_2;
	display z4_11_2;
	display x5_0_2;
	display z5_0_2;
	display x5_11_2;
	display z5_11_2;
	display x6_0_2;
	display z6_0_2;
	display x6_11_2;
	display z6_11_2;
	display x7_0_2;
	display z7_0_2;
	display x7_11_2;
	display z7_11_2;
	display x8_0_2;
	display z8_0_2;
	display x8_11_2;
	display z8_11_2;
	display x9_0_2;
	display z9_0_2;
	display x9_11_2;
	display z9_11_2;
	display x10_0_2;
	display z10_0_2;
	display x10_11_2;
	display z10_11_2;
	display x1_0_3;
	display z1_0_3;
	display x1_11_3;
	display z1_11_3;
	display x2_0_3;
	display z2_0_3;
	display x2_11_3;
	display z2_11_3;
	display x3_0_3;
	display z3_0_3;
	display x3_11_3;
	display z3_11_3;
	display x4_0_3;
	display z4_0_3;
	display x4_11_3;
	display z4_11_3;
	display x5_0_3;
	display z5_0_3;
	display x5_11_3;
	display z5_11_3;
	display x6_0_3;
	display z6_0_3;
	display x6_11_3;
	display z6_11_3;
	display x7_0_3;
	display z7_0_3;
	display x7_11_3;
	display z7_11_3;
	display x8_0_3;
	display z8_0_3;
	display x8_11_3;
	display z8_11_3;
	display x9_0_3;
	display z9_0_3;
	display x9_11_3;
	display z9_11_3;
	display x10_0_3;
	display z10_0_3;
	display x10_11_3;
	display z10_11_3;
	display x1_0_4;
	display z1_0_4;
	display x1_11_4;
	display z1_11_4;
	display x2_0_4;
	display z2_0_4;
	display x2_11_4;
	display z2_11_4;
	display x3_0_4;
	display z3_0_4;
	display x3_11_4;
	display z3_11_4;
	display x4_0_4;
	display z4_0_4;
	display x4_11_4;
	display z4_11_4;
	display x5_0_4;
	display z5_0_4;
	display x5_11_4;
	display z5_11_4;
	display x6_0_4;
	display z6_0_4;
	display x6_11_4;
	display z6_11_4;
	display x7_0_4;
	display z7_0_4;
	display x7_11_4;
	display z7_11_4;
	display x8_0_4;
	display z8_0_4;
	display x8_11_4;
	display z8_11_4;
	display x9_0_4;
	display z9_0_4;
	display x9_11_4;
	display z9_11_4;
	display x10_0_4;
	display z10_0_4;
	display x10_11_4;
	display z10_11_4;
	display x1_0_5;
	display z1_0_5;
	display x1_11_5;
	display z1_11_5;
	display x2_0_5;
	display z2_0_5;
	display x2_11_5;
	display z2_11_5;
	display x3_0_5;
	display z3_0_5;
	display x3_11_5;
	display z3_11_5;
	display x4_0_5;
	display z4_0_5;
	display x4_11_5;
	display z4_11_5;
	display x5_0_5;
	display z5_0_5;
	display x5_11_5;
	display z5_11_5;
	display x6_0_5;
	display z6_0_5;
	display x6_11_5;
	display z6_11_5;
	display x7_0_5;
	display z7_0_5;
	display x7_11_5;
	display z7_11_5;
	display x8_0_5;
	display z8_0_5;
	display x8_11_5;
	display z8_11_5;
	display x9_0_5;
	display z9_0_5;
	display x9_11_5;
	display z9_11_5;
	display x10_0_5;
	display z10_0_5;
	display x10_11_5;
	display z10_11_5;
	display x1_0_6;
	display z1_0_6;
	display x1_11_6;
	display z1_11_6;
	display x2_0_6;
	display z2_0_6;
	display x2_11_6;
	display z2_11_6;
	display x3_0_6;
	display z3_0_6;
	display x3_11_6;
	display z3_11_6;
	display x4_0_6;
	display z4_0_6;
	display x4_11_6;
	display z4_11_6;
	display x5_0_6;
	display z5_0_6;
	display x5_11_6;
	display z5_11_6;
	display x6_0_6;
	display z6_0_6;
	display x6_11_6;
	display z6_11_6;
	display x7_0_6;
	display z7_0_6;
	display x7_11_6;
	display z7_11_6;
	display x8_0_6;
	display z8_0_6;
	display x8_11_6;
	display z8_11_6;
	display x9_0_6;
	display z9_0_6;
	display x9_11_6;
	display z9_11_6;
	display x10_0_6;
	display z10_0_6;
	display x10_11_6;
	display z10_11_6;
	display x1_0_7;
	display z1_0_7;
	display x1_11_7;
	display z1_11_7;
	display x2_0_7;
	display z2_0_7;
	display x2_11_7;
	display z2_11_7;
	display x3_0_7;
	display z3_0_7;
	display x3_11_7;
	display z3_11_7;
	display x4_0_7;
	display z4_0_7;
	display x4_11_7;
	display z4_11_7;
	display x5_0_7;
	display z5_0_7;
	display x5_11_7;
	display z5_11_7;
	display x6_0_7;
	display z6_0_7;
	display x6_11_7;
	display z6_11_7;
	display x7_0_7;
	display z7_0_7;
	display x7_11_7;
	display z7_11_7;
	display x8_0_7;
	display z8_0_7;
	display x8_11_7;
	display z8_11_7;
	display x9_0_7;
	display z9_0_7;
	display x9_11_7;
	display z9_11_7;
	display x10_0_7;
	display z10_0_7;
	display x10_11_7;
	display z10_11_7;
	display x1_0_8;
	display z1_0_8;
	display x1_11_8;
	display z1_11_8;
	display x2_0_8;
	display z2_0_8;
	display x2_11_8;
	display z2_11_8;
	display x3_0_8;
	display z3_0_8;
	display x3_11_8;
	display z3_11_8;
	display x4_0_8;
	display z4_0_8;
	display x4_11_8;
	display z4_11_8;
	display x5_0_8;
	display z5_0_8;
	display x5_11_8;
	display z5_11_8;
	display x6_0_8;
	display z6_0_8;
	display x6_11_8;
	display z6_11_8;
	display x7_0_8;
	display z7_0_8;
	display x7_11_8;
	display z7_11_8;
	display x8_0_8;
	display z8_0_8;
	display x8_11_8;
	display z8_11_8;
	display x9_0_8;
	display z9_0_8;
	display x9_11_8;
	display z9_11_8;
	display x10_0_8;
	display z10_0_8;
	display x10_11_8;
	display z10_11_8;
	display x1_0_9;
	display z1_0_9;
	display x1_11_9;
	display z1_11_9;
	display x2_0_9;
	display z2_0_9;
	display x2_11_9;
	display z2_11_9;
	display x3_0_9;
	display z3_0_9;
	display x3_11_9;
	display z3_11_9;
	display x4_0_9;
	display z4_0_9;
	display x4_11_9;
	display z4_11_9;
	display x5_0_9;
	display z5_0_9;
	display x5_11_9;
	display z5_11_9;
	display x6_0_9;
	display z6_0_9;
	display x6_11_9;
	display z6_11_9;
	display x7_0_9;
	display z7_0_9;
	display x7_11_9;
	display z7_11_9;
	display x8_0_9;
	display z8_0_9;
	display x8_11_9;
	display z8_11_9;
	display x9_0_9;
	display z9_0_9;
	display x9_11_9;
	display z9_11_9;
	display x10_0_9;
	display z10_0_9;
	display x10_11_9;
	display z10_11_9;
	display x1_0_10;
	display z1_0_10;
	display x1_11_10;
	display z1_11_10;
	display x2_0_10;
	display z2_0_10;
	display x2_11_10;
	display z2_11_10;
	display x3_0_10;
	display z3_0_10;
	display x3_11_10;
	display z3_11_10;
	display x4_0_10;
	display z4_0_10;
	display x4_11_10;
	display z4_11_10;
	display x5_0_10;
	display z5_0_10;
	display x5_11_10;
	display z5_11_10;
	display x6_0_10;
	display z6_0_10;
	display x6_11_10;
	display z6_11_10;
	display x7_0_10;
	display z7_0_10;
	display x7_11_10;
	display z7_11_10;
	display x8_0_10;
	display z8_0_10;
	display x8_11_10;
	display z8_11_10;
	display x9_0_10;
	display z9_0_10;
	display x9_11_10;
	display z9_11_10;
	display x10_0_10;
	display z10_0_10;
	display x10_11_10;
	display z10_11_10;
	display x1_1_0;
	display y1_1_0;
	display x1_1_11;
	display y1_1_11;
	display x2_1_0;
	display y2_1_0;
	display x2_1_11;
	display y2_1_11;
	display x3_1_0;
	display y3_1_0;
	display x3_1_11;
	display y3_1_11;
	display x4_1_0;
	display y4_1_0;
	display x4_1_11;
	display y4_1_11;
	display x5_1_0;
	display y5_1_0;
	display x5_1_11;
	display y5_1_11;
	display x6_1_0;
	display y6_1_0;
	display x6_1_11;
	display y6_1_11;
	display x7_1_0;
	display y7_1_0;
	display x7_1_11;
	display y7_1_11;
	display x8_1_0;
	display y8_1_0;
	display x8_1_11;
	display y8_1_11;
	display x9_1_0;
	display y9_1_0;
	display x9_1_11;
	display y9_1_11;
	display x10_1_0;
	display y10_1_0;
	display x10_1_11;
	display y10_1_11;
	display x1_2_0;
	display y1_2_0;
	display x1_2_11;
	display y1_2_11;
	display x2_2_0;
	display y2_2_0;
	display x2_2_11;
	display y2_2_11;
	display x3_2_0;
	display y3_2_0;
	display x3_2_11;
	display y3_2_11;
	display x4_2_0;
	display y4_2_0;
	display x4_2_11;
	display y4_2_11;
	display x5_2_0;
	display y5_2_0;
	display x5_2_11;
	display y5_2_11;
	display x6_2_0;
	display y6_2_0;
	display x6_2_11;
	display y6_2_11;
	display x7_2_0;
	display y7_2_0;
	display x7_2_11;
	display y7_2_11;
	display x8_2_0;
	display y8_2_0;
	display x8_2_11;
	display y8_2_11;
	display x9_2_0;
	display y9_2_0;
	display x9_2_11;
	display y9_2_11;
	display x10_2_0;
	display y10_2_0;
	display x10_2_11;
	display y10_2_11;
	display x1_3_0;
	display y1_3_0;
	display x1_3_11;
	display y1_3_11;
	display x2_3_0;
	display y2_3_0;
	display x2_3_11;
	display y2_3_11;
	display x3_3_0;
	display y3_3_0;
	display x3_3_11;
	display y3_3_11;
	display x4_3_0;
	display y4_3_0;
	display x4_3_11;
	display y4_3_11;
	display x5_3_0;
	display y5_3_0;
	display x5_3_11;
	display y5_3_11;
	display x6_3_0;
	display y6_3_0;
	display x6_3_11;
	display y6_3_11;
	display x7_3_0;
	display y7_3_0;
	display x7_3_11;
	display y7_3_11;
	display x8_3_0;
	display y8_3_0;
	display x8_3_11;
	display y8_3_11;
	display x9_3_0;
	display y9_3_0;
	display x9_3_11;
	display y9_3_11;
	display x10_3_0;
	display y10_3_0;
	display x10_3_11;
	display y10_3_11;
	display x1_4_0;
	display y1_4_0;
	display x1_4_11;
	display y1_4_11;
	display x2_4_0;
	display y2_4_0;
	display x2_4_11;
	display y2_4_11;
	display x3_4_0;
	display y3_4_0;
	display x3_4_11;
	display y3_4_11;
	display x4_4_0;
	display y4_4_0;
	display x4_4_11;
	display y4_4_11;
	display x5_4_0;
	display y5_4_0;
	display x5_4_11;
	display y5_4_11;
	display x6_4_0;
	display y6_4_0;
	display x6_4_11;
	display y6_4_11;
	display x7_4_0;
	display y7_4_0;
	display x7_4_11;
	display y7_4_11;
	display x8_4_0;
	display y8_4_0;
	display x8_4_11;
	display y8_4_11;
	display x9_4_0;
	display y9_4_0;
	display x9_4_11;
	display y9_4_11;
	display x10_4_0;
	display y10_4_0;
	display x10_4_11;
	display y10_4_11;
	display x1_5_0;
	display y1_5_0;
	display x1_5_11;
	display y1_5_11;
	display x2_5_0;
	display y2_5_0;
	display x2_5_11;
	display y2_5_11;
	display x3_5_0;
	display y3_5_0;
	display x3_5_11;
	display y3_5_11;
	display x4_5_0;
	display y4_5_0;
	display x4_5_11;
	display y4_5_11;
	display x5_5_0;
	display y5_5_0;
	display x5_5_11;
	display y5_5_11;
	display x6_5_0;
	display y6_5_0;
	display x6_5_11;
	display y6_5_11;
	display x7_5_0;
	display y7_5_0;
	display x7_5_11;
	display y7_5_11;
	display x8_5_0;
	display y8_5_0;
	display x8_5_11;
	display y8_5_11;
	display x9_5_0;
	display y9_5_0;
	display x9_5_11;
	display y9_5_11;
	display x10_5_0;
	display y10_5_0;
	display x10_5_11;
	display y10_5_11;
	display x1_6_0;
	display y1_6_0;
	display x1_6_11;
	display y1_6_11;
	display x2_6_0;
	display y2_6_0;
	display x2_6_11;
	display y2_6_11;
	display x3_6_0;
	display y3_6_0;
	display x3_6_11;
	display y3_6_11;
	display x4_6_0;
	display y4_6_0;
	display x4_6_11;
	display y4_6_11;
	display x5_6_0;
	display y5_6_0;
	display x5_6_11;
	display y5_6_11;
	display x6_6_0;
	display y6_6_0;
	display x6_6_11;
	display y6_6_11;
	display x7_6_0;
	display y7_6_0;
	display x7_6_11;
	display y7_6_11;
	display x8_6_0;
	display y8_6_0;
	display x8_6_11;
	display y8_6_11;
	display x9_6_0;
	display y9_6_0;
	display x9_6_11;
	display y9_6_11;
	display x10_6_0;
	display y10_6_0;
	display x10_6_11;
	display y10_6_11;
	display x1_7_0;
	display y1_7_0;
	display x1_7_11;
	display y1_7_11;
	display x2_7_0;
	display y2_7_0;
	display x2_7_11;
	display y2_7_11;
	display x3_7_0;
	display y3_7_0;
	display x3_7_11;
	display y3_7_11;
	display x4_7_0;
	display y4_7_0;
	display x4_7_11;
	display y4_7_11;
	display x5_7_0;
	display y5_7_0;
	display x5_7_11;
	display y5_7_11;
	display x6_7_0;
	display y6_7_0;
	display x6_7_11;
	display y6_7_11;
	display x7_7_0;
	display y7_7_0;
	display x7_7_11;
	display y7_7_11;
	display x8_7_0;
	display y8_7_0;
	display x8_7_11;
	display y8_7_11;
	display x9_7_0;
	display y9_7_0;
	display x9_7_11;
	display y9_7_11;
	display x10_7_0;
	display y10_7_0;
	display x10_7_11;
	display y10_7_11;
	display x1_8_0;
	display y1_8_0;
	display x1_8_11;
	display y1_8_11;
	display x2_8_0;
	display y2_8_0;
	display x2_8_11;
	display y2_8_11;
	display x3_8_0;
	display y3_8_0;
	display x3_8_11;
	display y3_8_11;
	display x4_8_0;
	display y4_8_0;
	display x4_8_11;
	display y4_8_11;
	display x5_8_0;
	display y5_8_0;
	display x5_8_11;
	display y5_8_11;
	display x6_8_0;
	display y6_8_0;
	display x6_8_11;
	display y6_8_11;
	display x7_8_0;
	display y7_8_0;
	display x7_8_11;
	display y7_8_11;
	display x8_8_0;
	display y8_8_0;
	display x8_8_11;
	display y8_8_11;
	display x9_8_0;
	display y9_8_0;
	display x9_8_11;
	display y9_8_11;
	display x10_8_0;
	display y10_8_0;
	display x10_8_11;
	display y10_8_11;
	display x1_9_0;
	display y1_9_0;
	display x1_9_11;
	display y1_9_11;
	display x2_9_0;
	display y2_9_0;
	display x2_9_11;
	display y2_9_11;
	display x3_9_0;
	display y3_9_0;
	display x3_9_11;
	display y3_9_11;
	display x4_9_0;
	display y4_9_0;
	display x4_9_11;
	display y4_9_11;
	display x5_9_0;
	display y5_9_0;
	display x5_9_11;
	display y5_9_11;
	display x6_9_0;
	display y6_9_0;
	display x6_9_11;
	display y6_9_11;
	display x7_9_0;
	display y7_9_0;
	display x7_9_11;
	display y7_9_11;
	display x8_9_0;
	display y8_9_0;
	display x8_9_11;
	display y8_9_11;
	display x9_9_0;
	display y9_9_0;
	display x9_9_11;
	display y9_9_11;
	display x10_9_0;
	display y10_9_0;
	display x10_9_11;
	display y10_9_11;
	display x1_10_0;
	display y1_10_0;
	display x1_10_11;
	display y1_10_11;
	display x2_10_0;
	display y2_10_0;
	display x2_10_11;
	display y2_10_11;
	display x3_10_0;
	display y3_10_0;
	display x3_10_11;
	display y3_10_11;
	display x4_10_0;
	display y4_10_0;
	display x4_10_11;
	display y4_10_11;
	display x5_10_0;
	display y5_10_0;
	display x5_10_11;
	display y5_10_11;
	display x6_10_0;
	display y6_10_0;
	display x6_10_11;
	display y6_10_11;
	display x7_10_0;
	display y7_10_0;
	display x7_10_11;
	display y7_10_11;
	display x8_10_0;
	display y8_10_0;
	display x8_10_11;
	display y8_10_11;
	display x9_10_0;
	display y9_10_0;
	display x9_10_11;
	display y9_10_11;
	display x10_10_0;
	display y10_10_0;
	display x10_10_11;
	display y10_10_11;
display obj;
