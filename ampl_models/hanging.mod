#
#**************************
# SET UP THE INITIAL DATA *
#**************************
#   Problem :
#   *********
#   A catenary problem in 3 dimensions.  A rectangular grid is hung from its
#   4 corners under gravity.  The problem is to determine the resulting shape.
#   Source:  
#   an example in a talk by Nesterova and Vial, LLN, 1994.
#   SIF input: Ph. Toint, November 1994.
#   classification LQR2-AY-V-V
#   dimension of the grid
#IE NX                  3              $ n = 27
#IE NY                  3
#IE NX                  5              $ n = 90
#IE NY                  6
#IE NX                  20             $ n = 1800
#IE NY                  30
#   maximal X and Y distances
#   useful constants
#   Solution
#LO SOLTN(3,3)          -6.1184107487
#LO SOLTN(5,6)          -77.260229515
#LO SOLTN(10,10)        -620.17603242
param nx := 10;
param ny := 10;
param lx := 1.8;
param ly := 1.8;
param nxm1 := -1 + (10);
param nym1 := -1 + (10);
param lx2 := (1.8) * (1.8);
param ly2 := (1.8) * (1.8);
param rnx := 10.0;
param rny := 10.0;
param im1 := -1 + (10);
param rim1 := 9.0;
param jm1 := -1 + (10);
param rjm1 := 9.0;
param jp1 := 1 + (9);
param ip1 := 1 + (9);

var x1_1 >= 0.0 ,  <= 0.0;
var y1_1 >= 0.0 ,  <= 0.0;
var z1_1 >= 0.0 ,  <= 0.0;
var x1_2;
var y1_2 := 1.0;
var z1_2;
var x1_3;
var y1_3 := 2.0;
var z1_3;
var x1_4;
var y1_4 := 3.0;
var z1_4;
var x1_5;
var y1_5 := 4.0;
var z1_5;
var x1_6;
var y1_6 := 5.0;
var z1_6;
var x1_7;
var y1_7 := 6.0;
var z1_7;
var x1_8;
var y1_8 := 7.0;
var z1_8;
var x1_9;
var y1_9 := 8.0;
var z1_9;
var x1_10 >= 0.0 ,  <= 0.0;
var y1_10 >= 10.0 ,  <= 10.0 ,  := 9.0;
var z1_10 >= 0.0 ,  <= 0.0;
var x2_1 := 1.0;
var y2_1;
var z2_1;
var x2_2 := 1.0;
var y2_2 := 1.0;
var z2_2;
var x2_3 := 1.0;
var y2_3 := 2.0;
var z2_3;
var x2_4 := 1.0;
var y2_4 := 3.0;
var z2_4;
var x2_5 := 1.0;
var y2_5 := 4.0;
var z2_5;
var x2_6 := 1.0;
var y2_6 := 5.0;
var z2_6;
var x2_7 := 1.0;
var y2_7 := 6.0;
var z2_7;
var x2_8 := 1.0;
var y2_8 := 7.0;
var z2_8;
var x2_9 := 1.0;
var y2_9 := 8.0;
var z2_9;
var x2_10 := 1.0;
var y2_10 := 9.0;
var z2_10;
var x3_1 := 2.0;
var y3_1;
var z3_1;
var x3_2 := 2.0;
var y3_2 := 1.0;
var z3_2;
var x3_3 := 2.0;
var y3_3 := 2.0;
var z3_3;
var x3_4 := 2.0;
var y3_4 := 3.0;
var z3_4;
var x3_5 := 2.0;
var y3_5 := 4.0;
var z3_5;
var x3_6 := 2.0;
var y3_6 := 5.0;
var z3_6;
var x3_7 := 2.0;
var y3_7 := 6.0;
var z3_7;
var x3_8 := 2.0;
var y3_8 := 7.0;
var z3_8;
var x3_9 := 2.0;
var y3_9 := 8.0;
var z3_9;
var x3_10 := 2.0;
var y3_10 := 9.0;
var z3_10;
var x4_1 := 3.0;
var y4_1;
var z4_1;
var x4_2 := 3.0;
var y4_2 := 1.0;
var z4_2;
var x4_3 := 3.0;
var y4_3 := 2.0;
var z4_3;
var x4_4 := 3.0;
var y4_4 := 3.0;
var z4_4;
var x4_5 := 3.0;
var y4_5 := 4.0;
var z4_5;
var x4_6 := 3.0;
var y4_6 := 5.0;
var z4_6;
var x4_7 := 3.0;
var y4_7 := 6.0;
var z4_7;
var x4_8 := 3.0;
var y4_8 := 7.0;
var z4_8;
var x4_9 := 3.0;
var y4_9 := 8.0;
var z4_9;
var x4_10 := 3.0;
var y4_10 := 9.0;
var z4_10;
var x5_1 := 4.0;
var y5_1;
var z5_1;
var x5_2 := 4.0;
var y5_2 := 1.0;
var z5_2;
var x5_3 := 4.0;
var y5_3 := 2.0;
var z5_3;
var x5_4 := 4.0;
var y5_4 := 3.0;
var z5_4;
var x5_5 := 4.0;
var y5_5 := 4.0;
var z5_5;
var x5_6 := 4.0;
var y5_6 := 5.0;
var z5_6;
var x5_7 := 4.0;
var y5_7 := 6.0;
var z5_7;
var x5_8 := 4.0;
var y5_8 := 7.0;
var z5_8;
var x5_9 := 4.0;
var y5_9 := 8.0;
var z5_9;
var x5_10 := 4.0;
var y5_10 := 9.0;
var z5_10;
var x6_1 := 5.0;
var y6_1;
var z6_1;
var x6_2 := 5.0;
var y6_2 := 1.0;
var z6_2;
var x6_3 := 5.0;
var y6_3 := 2.0;
var z6_3;
var x6_4 := 5.0;
var y6_4 := 3.0;
var z6_4;
var x6_5 := 5.0;
var y6_5 := 4.0;
var z6_5;
var x6_6 := 5.0;
var y6_6 := 5.0;
var z6_6;
var x6_7 := 5.0;
var y6_7 := 6.0;
var z6_7;
var x6_8 := 5.0;
var y6_8 := 7.0;
var z6_8;
var x6_9 := 5.0;
var y6_9 := 8.0;
var z6_9;
var x6_10 := 5.0;
var y6_10 := 9.0;
var z6_10;
var x7_1 := 6.0;
var y7_1;
var z7_1;
var x7_2 := 6.0;
var y7_2 := 1.0;
var z7_2;
var x7_3 := 6.0;
var y7_3 := 2.0;
var z7_3;
var x7_4 := 6.0;
var y7_4 := 3.0;
var z7_4;
var x7_5 := 6.0;
var y7_5 := 4.0;
var z7_5;
var x7_6 := 6.0;
var y7_6 := 5.0;
var z7_6;
var x7_7 := 6.0;
var y7_7 := 6.0;
var z7_7;
var x7_8 := 6.0;
var y7_8 := 7.0;
var z7_8;
var x7_9 := 6.0;
var y7_9 := 8.0;
var z7_9;
var x7_10 := 6.0;
var y7_10 := 9.0;
var z7_10;
var x8_1 := 7.0;
var y8_1;
var z8_1;
var x8_2 := 7.0;
var y8_2 := 1.0;
var z8_2;
var x8_3 := 7.0;
var y8_3 := 2.0;
var z8_3;
var x8_4 := 7.0;
var y8_4 := 3.0;
var z8_4;
var x8_5 := 7.0;
var y8_5 := 4.0;
var z8_5;
var x8_6 := 7.0;
var y8_6 := 5.0;
var z8_6;
var x8_7 := 7.0;
var y8_7 := 6.0;
var z8_7;
var x8_8 := 7.0;
var y8_8 := 7.0;
var z8_8;
var x8_9 := 7.0;
var y8_9 := 8.0;
var z8_9;
var x8_10 := 7.0;
var y8_10 := 9.0;
var z8_10;
var x9_1 := 8.0;
var y9_1;
var z9_1;
var x9_2 := 8.0;
var y9_2 := 1.0;
var z9_2;
var x9_3 := 8.0;
var y9_3 := 2.0;
var z9_3;
var x9_4 := 8.0;
var y9_4 := 3.0;
var z9_4;
var x9_5 := 8.0;
var y9_5 := 4.0;
var z9_5;
var x9_6 := 8.0;
var y9_6 := 5.0;
var z9_6;
var x9_7 := 8.0;
var y9_7 := 6.0;
var z9_7;
var x9_8 := 8.0;
var y9_8 := 7.0;
var z9_8;
var x9_9 := 8.0;
var y9_9 := 8.0;
var z9_9;
var x9_10 := 8.0;
var y9_10 := 9.0;
var z9_10;
var x10_1 >= 10.0 ,  <= 10.0 ,  := 9.0;
var y10_1 >= 0.0 ,  <= 0.0;
var z10_1 >= 0.0 ,  <= 0.0;
var x10_2 := 9.0;
var y10_2 := 1.0;
var z10_2;
var x10_3 := 9.0;
var y10_3 := 2.0;
var z10_3;
var x10_4 := 9.0;
var y10_4 := 3.0;
var z10_4;
var x10_5 := 9.0;
var y10_5 := 4.0;
var z10_5;
var x10_6 := 9.0;
var y10_6 := 5.0;
var z10_6;
var x10_7 := 9.0;
var y10_7 := 6.0;
var z10_7;
var x10_8 := 9.0;
var y10_8 := 7.0;
var z10_8;
var x10_9 := 9.0;
var y10_9 := 8.0;
var z10_9;
var x10_10 >= 10.0 ,  <= 10.0 ,  := 9.0;
var y10_10 >= 10.0 ,  <= 10.0 ,  := 9.0;
var z10_10 >= 0.0 ,  <= 0.0;

minimize obj:
	z1_1 + z1_2 + z1_3 + z1_4 + z1_5 + z1_6 + z1_7 + z1_8 + z1_9 + z1_10 + z2_1 + z2_2 + z2_3 + z2_4 + z2_5 + z2_6 + z2_7 + z2_8 + z2_9 + z2_10 + z3_1 + z3_2 + z3_3 + z3_4 + z3_5 + z3_6 + z3_7 + z3_8 + z3_9 + z3_10 + z4_1 + z4_2 + z4_3 + z4_4 + z4_5 + z4_6 + z4_7 + z4_8 + z4_9 + z4_10 + z5_1 + z5_2 + z5_3 + z5_4 + z5_5 + z5_6 + z5_7 + z5_8 + z5_9 + z5_10 + z6_1 + z6_2 + z6_3 + z6_4 + z6_5 + z6_6 + z6_7 + z6_8 + z6_9 + z6_10 + z7_1 + z7_2 + z7_3 + z7_4 + z7_5 + z7_6 + z7_7 + z7_8 + z7_9 + z7_10 + z8_1 + z8_2 + z8_3 + z8_4 + z8_5 + z8_6 + z8_7 + z8_8 + z8_9 + z8_10 + z9_1 + z9_2 + z9_3 + z9_4 + z9_5 + z9_6 + z9_7 + z9_8 + z9_9 + z9_10 + z10_1 + z10_2 + z10_3 + z10_4 + z10_5 + z10_6 + z10_7 + z10_8 + z10_9 + z10_10;

subject to rc1_1:
	0 >= (x1_1-x1_2) * (x1_1-x1_2) + (y1_1-y1_2) * (y1_1-y1_2) + (z1_1-z1_2) * (z1_1-z1_2) - 3.24;
subject to rc1_2:
	0 >= (x1_2-x1_3) * (x1_2-x1_3) + (y1_2-y1_3) * (y1_2-y1_3) + (z1_2-z1_3) * (z1_2-z1_3) - 3.24;
subject to rc1_3:
	0 >= (x1_3-x1_4) * (x1_3-x1_4) + (y1_3-y1_4) * (y1_3-y1_4) + (z1_3-z1_4) * (z1_3-z1_4) - 3.24;
subject to rc1_4:
	0 >= (x1_4-x1_5) * (x1_4-x1_5) + (y1_4-y1_5) * (y1_4-y1_5) + (z1_4-z1_5) * (z1_4-z1_5) - 3.24;
subject to rc1_5:
	0 >= (x1_5-x1_6) * (x1_5-x1_6) + (y1_5-y1_6) * (y1_5-y1_6) + (z1_5-z1_6) * (z1_5-z1_6) - 3.24;
subject to rc1_6:
	0 >= (x1_6-x1_7) * (x1_6-x1_7) + (y1_6-y1_7) * (y1_6-y1_7) + (z1_6-z1_7) * (z1_6-z1_7) - 3.24;
subject to rc1_7:
	0 >= (x1_7-x1_8) * (x1_7-x1_8) + (y1_7-y1_8) * (y1_7-y1_8) + (z1_7-z1_8) * (z1_7-z1_8) - 3.24;
subject to rc1_8:
	0 >= (x1_8-x1_9) * (x1_8-x1_9) + (y1_8-y1_9) * (y1_8-y1_9) + (z1_8-z1_9) * (z1_8-z1_9) - 3.24;
subject to rc1_9:
	0 >= (x1_9-x1_10) * (x1_9-x1_10) + (y1_9-y1_10) * (y1_9-y1_10) + (z1_9-z1_10) * (z1_9-z1_10) - 3.24;
subject to rc2_1:
	0 >= (x2_1-x2_2) * (x2_1-x2_2) + (y2_1-y2_2) * (y2_1-y2_2) + (z2_1-z2_2) * (z2_1-z2_2) - 3.24;
subject to rc2_2:
	0 >= (x2_2-x2_3) * (x2_2-x2_3) + (y2_2-y2_3) * (y2_2-y2_3) + (z2_2-z2_3) * (z2_2-z2_3) - 3.24;
subject to rc2_3:
	0 >= (x2_3-x2_4) * (x2_3-x2_4) + (y2_3-y2_4) * (y2_3-y2_4) + (z2_3-z2_4) * (z2_3-z2_4) - 3.24;
subject to rc2_4:
	0 >= (x2_4-x2_5) * (x2_4-x2_5) + (y2_4-y2_5) * (y2_4-y2_5) + (z2_4-z2_5) * (z2_4-z2_5) - 3.24;
subject to rc2_5:
	0 >= (x2_5-x2_6) * (x2_5-x2_6) + (y2_5-y2_6) * (y2_5-y2_6) + (z2_5-z2_6) * (z2_5-z2_6) - 3.24;
subject to rc2_6:
	0 >= (x2_6-x2_7) * (x2_6-x2_7) + (y2_6-y2_7) * (y2_6-y2_7) + (z2_6-z2_7) * (z2_6-z2_7) - 3.24;
subject to rc2_7:
	0 >= (x2_7-x2_8) * (x2_7-x2_8) + (y2_7-y2_8) * (y2_7-y2_8) + (z2_7-z2_8) * (z2_7-z2_8) - 3.24;
subject to rc2_8:
	0 >= (x2_8-x2_9) * (x2_8-x2_9) + (y2_8-y2_9) * (y2_8-y2_9) + (z2_8-z2_9) * (z2_8-z2_9) - 3.24;
subject to rc2_9:
	0 >= (x2_9-x2_10) * (x2_9-x2_10) + (y2_9-y2_10) * (y2_9-y2_10) + (z2_9-z2_10) * (z2_9-z2_10) - 3.24;
subject to rc3_1:
	0 >= (x3_1-x3_2) * (x3_1-x3_2) + (y3_1-y3_2) * (y3_1-y3_2) + (z3_1-z3_2) * (z3_1-z3_2) - 3.24;
subject to rc3_2:
	0 >= (x3_2-x3_3) * (x3_2-x3_3) + (y3_2-y3_3) * (y3_2-y3_3) + (z3_2-z3_3) * (z3_2-z3_3) - 3.24;
subject to rc3_3:
	0 >= (x3_3-x3_4) * (x3_3-x3_4) + (y3_3-y3_4) * (y3_3-y3_4) + (z3_3-z3_4) * (z3_3-z3_4) - 3.24;
subject to rc3_4:
	0 >= (x3_4-x3_5) * (x3_4-x3_5) + (y3_4-y3_5) * (y3_4-y3_5) + (z3_4-z3_5) * (z3_4-z3_5) - 3.24;
subject to rc3_5:
	0 >= (x3_5-x3_6) * (x3_5-x3_6) + (y3_5-y3_6) * (y3_5-y3_6) + (z3_5-z3_6) * (z3_5-z3_6) - 3.24;
subject to rc3_6:
	0 >= (x3_6-x3_7) * (x3_6-x3_7) + (y3_6-y3_7) * (y3_6-y3_7) + (z3_6-z3_7) * (z3_6-z3_7) - 3.24;
subject to rc3_7:
	0 >= (x3_7-x3_8) * (x3_7-x3_8) + (y3_7-y3_8) * (y3_7-y3_8) + (z3_7-z3_8) * (z3_7-z3_8) - 3.24;
subject to rc3_8:
	0 >= (x3_8-x3_9) * (x3_8-x3_9) + (y3_8-y3_9) * (y3_8-y3_9) + (z3_8-z3_9) * (z3_8-z3_9) - 3.24;
subject to rc3_9:
	0 >= (x3_9-x3_10) * (x3_9-x3_10) + (y3_9-y3_10) * (y3_9-y3_10) + (z3_9-z3_10) * (z3_9-z3_10) - 3.24;
subject to rc4_1:
	0 >= (x4_1-x4_2) * (x4_1-x4_2) + (y4_1-y4_2) * (y4_1-y4_2) + (z4_1-z4_2) * (z4_1-z4_2) - 3.24;
subject to rc4_2:
	0 >= (x4_2-x4_3) * (x4_2-x4_3) + (y4_2-y4_3) * (y4_2-y4_3) + (z4_2-z4_3) * (z4_2-z4_3) - 3.24;
subject to rc4_3:
	0 >= (x4_3-x4_4) * (x4_3-x4_4) + (y4_3-y4_4) * (y4_3-y4_4) + (z4_3-z4_4) * (z4_3-z4_4) - 3.24;
subject to rc4_4:
	0 >= (x4_4-x4_5) * (x4_4-x4_5) + (y4_4-y4_5) * (y4_4-y4_5) + (z4_4-z4_5) * (z4_4-z4_5) - 3.24;
subject to rc4_5:
	0 >= (x4_5-x4_6) * (x4_5-x4_6) + (y4_5-y4_6) * (y4_5-y4_6) + (z4_5-z4_6) * (z4_5-z4_6) - 3.24;
subject to rc4_6:
	0 >= (x4_6-x4_7) * (x4_6-x4_7) + (y4_6-y4_7) * (y4_6-y4_7) + (z4_6-z4_7) * (z4_6-z4_7) - 3.24;
subject to rc4_7:
	0 >= (x4_7-x4_8) * (x4_7-x4_8) + (y4_7-y4_8) * (y4_7-y4_8) + (z4_7-z4_8) * (z4_7-z4_8) - 3.24;
subject to rc4_8:
	0 >= (x4_8-x4_9) * (x4_8-x4_9) + (y4_8-y4_9) * (y4_8-y4_9) + (z4_8-z4_9) * (z4_8-z4_9) - 3.24;
subject to rc4_9:
	0 >= (x4_9-x4_10) * (x4_9-x4_10) + (y4_9-y4_10) * (y4_9-y4_10) + (z4_9-z4_10) * (z4_9-z4_10) - 3.24;
subject to rc5_1:
	0 >= (x5_1-x5_2) * (x5_1-x5_2) + (y5_1-y5_2) * (y5_1-y5_2) + (z5_1-z5_2) * (z5_1-z5_2) - 3.24;
subject to rc5_2:
	0 >= (x5_2-x5_3) * (x5_2-x5_3) + (y5_2-y5_3) * (y5_2-y5_3) + (z5_2-z5_3) * (z5_2-z5_3) - 3.24;
subject to rc5_3:
	0 >= (x5_3-x5_4) * (x5_3-x5_4) + (y5_3-y5_4) * (y5_3-y5_4) + (z5_3-z5_4) * (z5_3-z5_4) - 3.24;
subject to rc5_4:
	0 >= (x5_4-x5_5) * (x5_4-x5_5) + (y5_4-y5_5) * (y5_4-y5_5) + (z5_4-z5_5) * (z5_4-z5_5) - 3.24;
subject to rc5_5:
	0 >= (x5_5-x5_6) * (x5_5-x5_6) + (y5_5-y5_6) * (y5_5-y5_6) + (z5_5-z5_6) * (z5_5-z5_6) - 3.24;
subject to rc5_6:
	0 >= (x5_6-x5_7) * (x5_6-x5_7) + (y5_6-y5_7) * (y5_6-y5_7) + (z5_6-z5_7) * (z5_6-z5_7) - 3.24;
subject to rc5_7:
	0 >= (x5_7-x5_8) * (x5_7-x5_8) + (y5_7-y5_8) * (y5_7-y5_8) + (z5_7-z5_8) * (z5_7-z5_8) - 3.24;
subject to rc5_8:
	0 >= (x5_8-x5_9) * (x5_8-x5_9) + (y5_8-y5_9) * (y5_8-y5_9) + (z5_8-z5_9) * (z5_8-z5_9) - 3.24;
subject to rc5_9:
	0 >= (x5_9-x5_10) * (x5_9-x5_10) + (y5_9-y5_10) * (y5_9-y5_10) + (z5_9-z5_10) * (z5_9-z5_10) - 3.24;
subject to rc6_1:
	0 >= (x6_1-x6_2) * (x6_1-x6_2) + (y6_1-y6_2) * (y6_1-y6_2) + (z6_1-z6_2) * (z6_1-z6_2) - 3.24;
subject to rc6_2:
	0 >= (x6_2-x6_3) * (x6_2-x6_3) + (y6_2-y6_3) * (y6_2-y6_3) + (z6_2-z6_3) * (z6_2-z6_3) - 3.24;
subject to rc6_3:
	0 >= (x6_3-x6_4) * (x6_3-x6_4) + (y6_3-y6_4) * (y6_3-y6_4) + (z6_3-z6_4) * (z6_3-z6_4) - 3.24;
subject to rc6_4:
	0 >= (x6_4-x6_5) * (x6_4-x6_5) + (y6_4-y6_5) * (y6_4-y6_5) + (z6_4-z6_5) * (z6_4-z6_5) - 3.24;
subject to rc6_5:
	0 >= (x6_5-x6_6) * (x6_5-x6_6) + (y6_5-y6_6) * (y6_5-y6_6) + (z6_5-z6_6) * (z6_5-z6_6) - 3.24;
subject to rc6_6:
	0 >= (x6_6-x6_7) * (x6_6-x6_7) + (y6_6-y6_7) * (y6_6-y6_7) + (z6_6-z6_7) * (z6_6-z6_7) - 3.24;
subject to rc6_7:
	0 >= (x6_7-x6_8) * (x6_7-x6_8) + (y6_7-y6_8) * (y6_7-y6_8) + (z6_7-z6_8) * (z6_7-z6_8) - 3.24;
subject to rc6_8:
	0 >= (x6_8-x6_9) * (x6_8-x6_9) + (y6_8-y6_9) * (y6_8-y6_9) + (z6_8-z6_9) * (z6_8-z6_9) - 3.24;
subject to rc6_9:
	0 >= (x6_9-x6_10) * (x6_9-x6_10) + (y6_9-y6_10) * (y6_9-y6_10) + (z6_9-z6_10) * (z6_9-z6_10) - 3.24;
subject to rc7_1:
	0 >= (x7_1-x7_2) * (x7_1-x7_2) + (y7_1-y7_2) * (y7_1-y7_2) + (z7_1-z7_2) * (z7_1-z7_2) - 3.24;
subject to rc7_2:
	0 >= (x7_2-x7_3) * (x7_2-x7_3) + (y7_2-y7_3) * (y7_2-y7_3) + (z7_2-z7_3) * (z7_2-z7_3) - 3.24;
subject to rc7_3:
	0 >= (x7_3-x7_4) * (x7_3-x7_4) + (y7_3-y7_4) * (y7_3-y7_4) + (z7_3-z7_4) * (z7_3-z7_4) - 3.24;
subject to rc7_4:
	0 >= (x7_4-x7_5) * (x7_4-x7_5) + (y7_4-y7_5) * (y7_4-y7_5) + (z7_4-z7_5) * (z7_4-z7_5) - 3.24;
subject to rc7_5:
	0 >= (x7_5-x7_6) * (x7_5-x7_6) + (y7_5-y7_6) * (y7_5-y7_6) + (z7_5-z7_6) * (z7_5-z7_6) - 3.24;
subject to rc7_6:
	0 >= (x7_6-x7_7) * (x7_6-x7_7) + (y7_6-y7_7) * (y7_6-y7_7) + (z7_6-z7_7) * (z7_6-z7_7) - 3.24;
subject to rc7_7:
	0 >= (x7_7-x7_8) * (x7_7-x7_8) + (y7_7-y7_8) * (y7_7-y7_8) + (z7_7-z7_8) * (z7_7-z7_8) - 3.24;
subject to rc7_8:
	0 >= (x7_8-x7_9) * (x7_8-x7_9) + (y7_8-y7_9) * (y7_8-y7_9) + (z7_8-z7_9) * (z7_8-z7_9) - 3.24;
subject to rc7_9:
	0 >= (x7_9-x7_10) * (x7_9-x7_10) + (y7_9-y7_10) * (y7_9-y7_10) + (z7_9-z7_10) * (z7_9-z7_10) - 3.24;
subject to rc8_1:
	0 >= (x8_1-x8_2) * (x8_1-x8_2) + (y8_1-y8_2) * (y8_1-y8_2) + (z8_1-z8_2) * (z8_1-z8_2) - 3.24;
subject to rc8_2:
	0 >= (x8_2-x8_3) * (x8_2-x8_3) + (y8_2-y8_3) * (y8_2-y8_3) + (z8_2-z8_3) * (z8_2-z8_3) - 3.24;
subject to rc8_3:
	0 >= (x8_3-x8_4) * (x8_3-x8_4) + (y8_3-y8_4) * (y8_3-y8_4) + (z8_3-z8_4) * (z8_3-z8_4) - 3.24;
subject to rc8_4:
	0 >= (x8_4-x8_5) * (x8_4-x8_5) + (y8_4-y8_5) * (y8_4-y8_5) + (z8_4-z8_5) * (z8_4-z8_5) - 3.24;
subject to rc8_5:
	0 >= (x8_5-x8_6) * (x8_5-x8_6) + (y8_5-y8_6) * (y8_5-y8_6) + (z8_5-z8_6) * (z8_5-z8_6) - 3.24;
subject to rc8_6:
	0 >= (x8_6-x8_7) * (x8_6-x8_7) + (y8_6-y8_7) * (y8_6-y8_7) + (z8_6-z8_7) * (z8_6-z8_7) - 3.24;
subject to rc8_7:
	0 >= (x8_7-x8_8) * (x8_7-x8_8) + (y8_7-y8_8) * (y8_7-y8_8) + (z8_7-z8_8) * (z8_7-z8_8) - 3.24;
subject to rc8_8:
	0 >= (x8_8-x8_9) * (x8_8-x8_9) + (y8_8-y8_9) * (y8_8-y8_9) + (z8_8-z8_9) * (z8_8-z8_9) - 3.24;
subject to rc8_9:
	0 >= (x8_9-x8_10) * (x8_9-x8_10) + (y8_9-y8_10) * (y8_9-y8_10) + (z8_9-z8_10) * (z8_9-z8_10) - 3.24;
subject to rc9_1:
	0 >= (x9_1-x9_2) * (x9_1-x9_2) + (y9_1-y9_2) * (y9_1-y9_2) + (z9_1-z9_2) * (z9_1-z9_2) - 3.24;
subject to rc9_2:
	0 >= (x9_2-x9_3) * (x9_2-x9_3) + (y9_2-y9_3) * (y9_2-y9_3) + (z9_2-z9_3) * (z9_2-z9_3) - 3.24;
subject to rc9_3:
	0 >= (x9_3-x9_4) * (x9_3-x9_4) + (y9_3-y9_4) * (y9_3-y9_4) + (z9_3-z9_4) * (z9_3-z9_4) - 3.24;
subject to rc9_4:
	0 >= (x9_4-x9_5) * (x9_4-x9_5) + (y9_4-y9_5) * (y9_4-y9_5) + (z9_4-z9_5) * (z9_4-z9_5) - 3.24;
subject to rc9_5:
	0 >= (x9_5-x9_6) * (x9_5-x9_6) + (y9_5-y9_6) * (y9_5-y9_6) + (z9_5-z9_6) * (z9_5-z9_6) - 3.24;
subject to rc9_6:
	0 >= (x9_6-x9_7) * (x9_6-x9_7) + (y9_6-y9_7) * (y9_6-y9_7) + (z9_6-z9_7) * (z9_6-z9_7) - 3.24;
subject to rc9_7:
	0 >= (x9_7-x9_8) * (x9_7-x9_8) + (y9_7-y9_8) * (y9_7-y9_8) + (z9_7-z9_8) * (z9_7-z9_8) - 3.24;
subject to rc9_8:
	0 >= (x9_8-x9_9) * (x9_8-x9_9) + (y9_8-y9_9) * (y9_8-y9_9) + (z9_8-z9_9) * (z9_8-z9_9) - 3.24;
subject to rc9_9:
	0 >= (x9_9-x9_10) * (x9_9-x9_10) + (y9_9-y9_10) * (y9_9-y9_10) + (z9_9-z9_10) * (z9_9-z9_10) - 3.24;
subject to rc10_1:
	0 >= (x10_1-x10_2) * (x10_1-x10_2) + (y10_1-y10_2) * (y10_1-y10_2) + (z10_1-z10_2) * (z10_1-z10_2) - 3.24;
subject to rc10_2:
	0 >= (x10_2-x10_3) * (x10_2-x10_3) + (y10_2-y10_3) * (y10_2-y10_3) + (z10_2-z10_3) * (z10_2-z10_3) - 3.24;
subject to rc10_3:
	0 >= (x10_3-x10_4) * (x10_3-x10_4) + (y10_3-y10_4) * (y10_3-y10_4) + (z10_3-z10_4) * (z10_3-z10_4) - 3.24;
subject to rc10_4:
	0 >= (x10_4-x10_5) * (x10_4-x10_5) + (y10_4-y10_5) * (y10_4-y10_5) + (z10_4-z10_5) * (z10_4-z10_5) - 3.24;
subject to rc10_5:
	0 >= (x10_5-x10_6) * (x10_5-x10_6) + (y10_5-y10_6) * (y10_5-y10_6) + (z10_5-z10_6) * (z10_5-z10_6) - 3.24;
subject to rc10_6:
	0 >= (x10_6-x10_7) * (x10_6-x10_7) + (y10_6-y10_7) * (y10_6-y10_7) + (z10_6-z10_7) * (z10_6-z10_7) - 3.24;
subject to rc10_7:
	0 >= (x10_7-x10_8) * (x10_7-x10_8) + (y10_7-y10_8) * (y10_7-y10_8) + (z10_7-z10_8) * (z10_7-z10_8) - 3.24;
subject to rc10_8:
	0 >= (x10_8-x10_9) * (x10_8-x10_9) + (y10_8-y10_9) * (y10_8-y10_9) + (z10_8-z10_9) * (z10_8-z10_9) - 3.24;
subject to rc10_9:
	0 >= (x10_9-x10_10) * (x10_9-x10_10) + (y10_9-y10_10) * (y10_9-y10_10) + (z10_9-z10_10) * (z10_9-z10_10) - 3.24;
subject to dc1_1:
	0 >= (x1_1-x2_1) * (x1_1-x2_1) + (y1_1-y2_1) * (y1_1-y2_1) + (z1_1-z2_1) * (z1_1-z2_1) - 3.24;
subject to dc1_2:
	0 >= (x1_2-x2_2) * (x1_2-x2_2) + (y1_2-y2_2) * (y1_2-y2_2) + (z1_2-z2_2) * (z1_2-z2_2) - 3.24;
subject to dc1_3:
	0 >= (x1_3-x2_3) * (x1_3-x2_3) + (y1_3-y2_3) * (y1_3-y2_3) + (z1_3-z2_3) * (z1_3-z2_3) - 3.24;
subject to dc1_4:
	0 >= (x1_4-x2_4) * (x1_4-x2_4) + (y1_4-y2_4) * (y1_4-y2_4) + (z1_4-z2_4) * (z1_4-z2_4) - 3.24;
subject to dc1_5:
	0 >= (x1_5-x2_5) * (x1_5-x2_5) + (y1_5-y2_5) * (y1_5-y2_5) + (z1_5-z2_5) * (z1_5-z2_5) - 3.24;
subject to dc1_6:
	0 >= (x1_6-x2_6) * (x1_6-x2_6) + (y1_6-y2_6) * (y1_6-y2_6) + (z1_6-z2_6) * (z1_6-z2_6) - 3.24;
subject to dc1_7:
	0 >= (x1_7-x2_7) * (x1_7-x2_7) + (y1_7-y2_7) * (y1_7-y2_7) + (z1_7-z2_7) * (z1_7-z2_7) - 3.24;
subject to dc1_8:
	0 >= (x1_8-x2_8) * (x1_8-x2_8) + (y1_8-y2_8) * (y1_8-y2_8) + (z1_8-z2_8) * (z1_8-z2_8) - 3.24;
subject to dc1_9:
	0 >= (x1_9-x2_9) * (x1_9-x2_9) + (y1_9-y2_9) * (y1_9-y2_9) + (z1_9-z2_9) * (z1_9-z2_9) - 3.24;
subject to dc1_10:
	0 >= (x1_10-x2_10) * (x1_10-x2_10) + (y1_10-y2_10) * (y1_10-y2_10) + (z1_10-z2_10) * (z1_10-z2_10) - 3.24;
subject to dc2_1:
	0 >= (x2_1-x3_1) * (x2_1-x3_1) + (y2_1-y3_1) * (y2_1-y3_1) + (z2_1-z3_1) * (z2_1-z3_1) - 3.24;
subject to dc2_2:
	0 >= (x2_2-x3_2) * (x2_2-x3_2) + (y2_2-y3_2) * (y2_2-y3_2) + (z2_2-z3_2) * (z2_2-z3_2) - 3.24;
subject to dc2_3:
	0 >= (x2_3-x3_3) * (x2_3-x3_3) + (y2_3-y3_3) * (y2_3-y3_3) + (z2_3-z3_3) * (z2_3-z3_3) - 3.24;
subject to dc2_4:
	0 >= (x2_4-x3_4) * (x2_4-x3_4) + (y2_4-y3_4) * (y2_4-y3_4) + (z2_4-z3_4) * (z2_4-z3_4) - 3.24;
subject to dc2_5:
	0 >= (x2_5-x3_5) * (x2_5-x3_5) + (y2_5-y3_5) * (y2_5-y3_5) + (z2_5-z3_5) * (z2_5-z3_5) - 3.24;
subject to dc2_6:
	0 >= (x2_6-x3_6) * (x2_6-x3_6) + (y2_6-y3_6) * (y2_6-y3_6) + (z2_6-z3_6) * (z2_6-z3_6) - 3.24;
subject to dc2_7:
	0 >= (x2_7-x3_7) * (x2_7-x3_7) + (y2_7-y3_7) * (y2_7-y3_7) + (z2_7-z3_7) * (z2_7-z3_7) - 3.24;
subject to dc2_8:
	0 >= (x2_8-x3_8) * (x2_8-x3_8) + (y2_8-y3_8) * (y2_8-y3_8) + (z2_8-z3_8) * (z2_8-z3_8) - 3.24;
subject to dc2_9:
	0 >= (x2_9-x3_9) * (x2_9-x3_9) + (y2_9-y3_9) * (y2_9-y3_9) + (z2_9-z3_9) * (z2_9-z3_9) - 3.24;
subject to dc2_10:
	0 >= (x2_10-x3_10) * (x2_10-x3_10) + (y2_10-y3_10) * (y2_10-y3_10) + (z2_10-z3_10) * (z2_10-z3_10) - 3.24;
subject to dc3_1:
	0 >= (x3_1-x4_1) * (x3_1-x4_1) + (y3_1-y4_1) * (y3_1-y4_1) + (z3_1-z4_1) * (z3_1-z4_1) - 3.24;
subject to dc3_2:
	0 >= (x3_2-x4_2) * (x3_2-x4_2) + (y3_2-y4_2) * (y3_2-y4_2) + (z3_2-z4_2) * (z3_2-z4_2) - 3.24;
subject to dc3_3:
	0 >= (x3_3-x4_3) * (x3_3-x4_3) + (y3_3-y4_3) * (y3_3-y4_3) + (z3_3-z4_3) * (z3_3-z4_3) - 3.24;
subject to dc3_4:
	0 >= (x3_4-x4_4) * (x3_4-x4_4) + (y3_4-y4_4) * (y3_4-y4_4) + (z3_4-z4_4) * (z3_4-z4_4) - 3.24;
subject to dc3_5:
	0 >= (x3_5-x4_5) * (x3_5-x4_5) + (y3_5-y4_5) * (y3_5-y4_5) + (z3_5-z4_5) * (z3_5-z4_5) - 3.24;
subject to dc3_6:
	0 >= (x3_6-x4_6) * (x3_6-x4_6) + (y3_6-y4_6) * (y3_6-y4_6) + (z3_6-z4_6) * (z3_6-z4_6) - 3.24;
subject to dc3_7:
	0 >= (x3_7-x4_7) * (x3_7-x4_7) + (y3_7-y4_7) * (y3_7-y4_7) + (z3_7-z4_7) * (z3_7-z4_7) - 3.24;
subject to dc3_8:
	0 >= (x3_8-x4_8) * (x3_8-x4_8) + (y3_8-y4_8) * (y3_8-y4_8) + (z3_8-z4_8) * (z3_8-z4_8) - 3.24;
subject to dc3_9:
	0 >= (x3_9-x4_9) * (x3_9-x4_9) + (y3_9-y4_9) * (y3_9-y4_9) + (z3_9-z4_9) * (z3_9-z4_9) - 3.24;
subject to dc3_10:
	0 >= (x3_10-x4_10) * (x3_10-x4_10) + (y3_10-y4_10) * (y3_10-y4_10) + (z3_10-z4_10) * (z3_10-z4_10) - 3.24;
subject to dc4_1:
	0 >= (x4_1-x5_1) * (x4_1-x5_1) + (y4_1-y5_1) * (y4_1-y5_1) + (z4_1-z5_1) * (z4_1-z5_1) - 3.24;
subject to dc4_2:
	0 >= (x4_2-x5_2) * (x4_2-x5_2) + (y4_2-y5_2) * (y4_2-y5_2) + (z4_2-z5_2) * (z4_2-z5_2) - 3.24;
subject to dc4_3:
	0 >= (x4_3-x5_3) * (x4_3-x5_3) + (y4_3-y5_3) * (y4_3-y5_3) + (z4_3-z5_3) * (z4_3-z5_3) - 3.24;
subject to dc4_4:
	0 >= (x4_4-x5_4) * (x4_4-x5_4) + (y4_4-y5_4) * (y4_4-y5_4) + (z4_4-z5_4) * (z4_4-z5_4) - 3.24;
subject to dc4_5:
	0 >= (x4_5-x5_5) * (x4_5-x5_5) + (y4_5-y5_5) * (y4_5-y5_5) + (z4_5-z5_5) * (z4_5-z5_5) - 3.24;
subject to dc4_6:
	0 >= (x4_6-x5_6) * (x4_6-x5_6) + (y4_6-y5_6) * (y4_6-y5_6) + (z4_6-z5_6) * (z4_6-z5_6) - 3.24;
subject to dc4_7:
	0 >= (x4_7-x5_7) * (x4_7-x5_7) + (y4_7-y5_7) * (y4_7-y5_7) + (z4_7-z5_7) * (z4_7-z5_7) - 3.24;
subject to dc4_8:
	0 >= (x4_8-x5_8) * (x4_8-x5_8) + (y4_8-y5_8) * (y4_8-y5_8) + (z4_8-z5_8) * (z4_8-z5_8) - 3.24;
subject to dc4_9:
	0 >= (x4_9-x5_9) * (x4_9-x5_9) + (y4_9-y5_9) * (y4_9-y5_9) + (z4_9-z5_9) * (z4_9-z5_9) - 3.24;
subject to dc4_10:
	0 >= (x4_10-x5_10) * (x4_10-x5_10) + (y4_10-y5_10) * (y4_10-y5_10) + (z4_10-z5_10) * (z4_10-z5_10) - 3.24;
subject to dc5_1:
	0 >= (x5_1-x6_1) * (x5_1-x6_1) + (y5_1-y6_1) * (y5_1-y6_1) + (z5_1-z6_1) * (z5_1-z6_1) - 3.24;
subject to dc5_2:
	0 >= (x5_2-x6_2) * (x5_2-x6_2) + (y5_2-y6_2) * (y5_2-y6_2) + (z5_2-z6_2) * (z5_2-z6_2) - 3.24;
subject to dc5_3:
	0 >= (x5_3-x6_3) * (x5_3-x6_3) + (y5_3-y6_3) * (y5_3-y6_3) + (z5_3-z6_3) * (z5_3-z6_3) - 3.24;
subject to dc5_4:
	0 >= (x5_4-x6_4) * (x5_4-x6_4) + (y5_4-y6_4) * (y5_4-y6_4) + (z5_4-z6_4) * (z5_4-z6_4) - 3.24;
subject to dc5_5:
	0 >= (x5_5-x6_5) * (x5_5-x6_5) + (y5_5-y6_5) * (y5_5-y6_5) + (z5_5-z6_5) * (z5_5-z6_5) - 3.24;
subject to dc5_6:
	0 >= (x5_6-x6_6) * (x5_6-x6_6) + (y5_6-y6_6) * (y5_6-y6_6) + (z5_6-z6_6) * (z5_6-z6_6) - 3.24;
subject to dc5_7:
	0 >= (x5_7-x6_7) * (x5_7-x6_7) + (y5_7-y6_7) * (y5_7-y6_7) + (z5_7-z6_7) * (z5_7-z6_7) - 3.24;
subject to dc5_8:
	0 >= (x5_8-x6_8) * (x5_8-x6_8) + (y5_8-y6_8) * (y5_8-y6_8) + (z5_8-z6_8) * (z5_8-z6_8) - 3.24;
subject to dc5_9:
	0 >= (x5_9-x6_9) * (x5_9-x6_9) + (y5_9-y6_9) * (y5_9-y6_9) + (z5_9-z6_9) * (z5_9-z6_9) - 3.24;
subject to dc5_10:
	0 >= (x5_10-x6_10) * (x5_10-x6_10) + (y5_10-y6_10) * (y5_10-y6_10) + (z5_10-z6_10) * (z5_10-z6_10) - 3.24;
subject to dc6_1:
	0 >= (x6_1-x7_1) * (x6_1-x7_1) + (y6_1-y7_1) * (y6_1-y7_1) + (z6_1-z7_1) * (z6_1-z7_1) - 3.24;
subject to dc6_2:
	0 >= (x6_2-x7_2) * (x6_2-x7_2) + (y6_2-y7_2) * (y6_2-y7_2) + (z6_2-z7_2) * (z6_2-z7_2) - 3.24;
subject to dc6_3:
	0 >= (x6_3-x7_3) * (x6_3-x7_3) + (y6_3-y7_3) * (y6_3-y7_3) + (z6_3-z7_3) * (z6_3-z7_3) - 3.24;
subject to dc6_4:
	0 >= (x6_4-x7_4) * (x6_4-x7_4) + (y6_4-y7_4) * (y6_4-y7_4) + (z6_4-z7_4) * (z6_4-z7_4) - 3.24;
subject to dc6_5:
	0 >= (x6_5-x7_5) * (x6_5-x7_5) + (y6_5-y7_5) * (y6_5-y7_5) + (z6_5-z7_5) * (z6_5-z7_5) - 3.24;
subject to dc6_6:
	0 >= (x6_6-x7_6) * (x6_6-x7_6) + (y6_6-y7_6) * (y6_6-y7_6) + (z6_6-z7_6) * (z6_6-z7_6) - 3.24;
subject to dc6_7:
	0 >= (x6_7-x7_7) * (x6_7-x7_7) + (y6_7-y7_7) * (y6_7-y7_7) + (z6_7-z7_7) * (z6_7-z7_7) - 3.24;
subject to dc6_8:
	0 >= (x6_8-x7_8) * (x6_8-x7_8) + (y6_8-y7_8) * (y6_8-y7_8) + (z6_8-z7_8) * (z6_8-z7_8) - 3.24;
subject to dc6_9:
	0 >= (x6_9-x7_9) * (x6_9-x7_9) + (y6_9-y7_9) * (y6_9-y7_9) + (z6_9-z7_9) * (z6_9-z7_9) - 3.24;
subject to dc6_10:
	0 >= (x6_10-x7_10) * (x6_10-x7_10) + (y6_10-y7_10) * (y6_10-y7_10) + (z6_10-z7_10) * (z6_10-z7_10) - 3.24;
subject to dc7_1:
	0 >= (x7_1-x8_1) * (x7_1-x8_1) + (y7_1-y8_1) * (y7_1-y8_1) + (z7_1-z8_1) * (z7_1-z8_1) - 3.24;
subject to dc7_2:
	0 >= (x7_2-x8_2) * (x7_2-x8_2) + (y7_2-y8_2) * (y7_2-y8_2) + (z7_2-z8_2) * (z7_2-z8_2) - 3.24;
subject to dc7_3:
	0 >= (x7_3-x8_3) * (x7_3-x8_3) + (y7_3-y8_3) * (y7_3-y8_3) + (z7_3-z8_3) * (z7_3-z8_3) - 3.24;
subject to dc7_4:
	0 >= (x7_4-x8_4) * (x7_4-x8_4) + (y7_4-y8_4) * (y7_4-y8_4) + (z7_4-z8_4) * (z7_4-z8_4) - 3.24;
subject to dc7_5:
	0 >= (x7_5-x8_5) * (x7_5-x8_5) + (y7_5-y8_5) * (y7_5-y8_5) + (z7_5-z8_5) * (z7_5-z8_5) - 3.24;
subject to dc7_6:
	0 >= (x7_6-x8_6) * (x7_6-x8_6) + (y7_6-y8_6) * (y7_6-y8_6) + (z7_6-z8_6) * (z7_6-z8_6) - 3.24;
subject to dc7_7:
	0 >= (x7_7-x8_7) * (x7_7-x8_7) + (y7_7-y8_7) * (y7_7-y8_7) + (z7_7-z8_7) * (z7_7-z8_7) - 3.24;
subject to dc7_8:
	0 >= (x7_8-x8_8) * (x7_8-x8_8) + (y7_8-y8_8) * (y7_8-y8_8) + (z7_8-z8_8) * (z7_8-z8_8) - 3.24;
subject to dc7_9:
	0 >= (x7_9-x8_9) * (x7_9-x8_9) + (y7_9-y8_9) * (y7_9-y8_9) + (z7_9-z8_9) * (z7_9-z8_9) - 3.24;
subject to dc7_10:
	0 >= (x7_10-x8_10) * (x7_10-x8_10) + (y7_10-y8_10) * (y7_10-y8_10) + (z7_10-z8_10) * (z7_10-z8_10) - 3.24;
subject to dc8_1:
	0 >= (x8_1-x9_1) * (x8_1-x9_1) + (y8_1-y9_1) * (y8_1-y9_1) + (z8_1-z9_1) * (z8_1-z9_1) - 3.24;
subject to dc8_2:
	0 >= (x8_2-x9_2) * (x8_2-x9_2) + (y8_2-y9_2) * (y8_2-y9_2) + (z8_2-z9_2) * (z8_2-z9_2) - 3.24;
subject to dc8_3:
	0 >= (x8_3-x9_3) * (x8_3-x9_3) + (y8_3-y9_3) * (y8_3-y9_3) + (z8_3-z9_3) * (z8_3-z9_3) - 3.24;
subject to dc8_4:
	0 >= (x8_4-x9_4) * (x8_4-x9_4) + (y8_4-y9_4) * (y8_4-y9_4) + (z8_4-z9_4) * (z8_4-z9_4) - 3.24;
subject to dc8_5:
	0 >= (x8_5-x9_5) * (x8_5-x9_5) + (y8_5-y9_5) * (y8_5-y9_5) + (z8_5-z9_5) * (z8_5-z9_5) - 3.24;
subject to dc8_6:
	0 >= (x8_6-x9_6) * (x8_6-x9_6) + (y8_6-y9_6) * (y8_6-y9_6) + (z8_6-z9_6) * (z8_6-z9_6) - 3.24;
subject to dc8_7:
	0 >= (x8_7-x9_7) * (x8_7-x9_7) + (y8_7-y9_7) * (y8_7-y9_7) + (z8_7-z9_7) * (z8_7-z9_7) - 3.24;
subject to dc8_8:
	0 >= (x8_8-x9_8) * (x8_8-x9_8) + (y8_8-y9_8) * (y8_8-y9_8) + (z8_8-z9_8) * (z8_8-z9_8) - 3.24;
subject to dc8_9:
	0 >= (x8_9-x9_9) * (x8_9-x9_9) + (y8_9-y9_9) * (y8_9-y9_9) + (z8_9-z9_9) * (z8_9-z9_9) - 3.24;
subject to dc8_10:
	0 >= (x8_10-x9_10) * (x8_10-x9_10) + (y8_10-y9_10) * (y8_10-y9_10) + (z8_10-z9_10) * (z8_10-z9_10) - 3.24;
subject to dc9_1:
	0 >= (x9_1-x10_1) * (x9_1-x10_1) + (y9_1-y10_1) * (y9_1-y10_1) + (z9_1-z10_1) * (z9_1-z10_1) - 3.24;
subject to dc9_2:
	0 >= (x9_2-x10_2) * (x9_2-x10_2) + (y9_2-y10_2) * (y9_2-y10_2) + (z9_2-z10_2) * (z9_2-z10_2) - 3.24;
subject to dc9_3:
	0 >= (x9_3-x10_3) * (x9_3-x10_3) + (y9_3-y10_3) * (y9_3-y10_3) + (z9_3-z10_3) * (z9_3-z10_3) - 3.24;
subject to dc9_4:
	0 >= (x9_4-x10_4) * (x9_4-x10_4) + (y9_4-y10_4) * (y9_4-y10_4) + (z9_4-z10_4) * (z9_4-z10_4) - 3.24;
subject to dc9_5:
	0 >= (x9_5-x10_5) * (x9_5-x10_5) + (y9_5-y10_5) * (y9_5-y10_5) + (z9_5-z10_5) * (z9_5-z10_5) - 3.24;
subject to dc9_6:
	0 >= (x9_6-x10_6) * (x9_6-x10_6) + (y9_6-y10_6) * (y9_6-y10_6) + (z9_6-z10_6) * (z9_6-z10_6) - 3.24;
subject to dc9_7:
	0 >= (x9_7-x10_7) * (x9_7-x10_7) + (y9_7-y10_7) * (y9_7-y10_7) + (z9_7-z10_7) * (z9_7-z10_7) - 3.24;
subject to dc9_8:
	0 >= (x9_8-x10_8) * (x9_8-x10_8) + (y9_8-y10_8) * (y9_8-y10_8) + (z9_8-z10_8) * (z9_8-z10_8) - 3.24;
subject to dc9_9:
	0 >= (x9_9-x10_9) * (x9_9-x10_9) + (y9_9-y10_9) * (y9_9-y10_9) + (z9_9-z10_9) * (z9_9-z10_9) - 3.24;
subject to dc9_10:
	0 >= (x9_10-x10_10) * (x9_10-x10_10) + (y9_10-y10_10) * (y9_10-y10_10) + (z9_10-z10_10) * (z9_10-z10_10) - 3.24;

solve;
display x1_1;
display y1_1;
display z1_1;
display x1_2;
display y1_2;
display z1_2;
display x1_3;
display y1_3;
display z1_3;
display x1_4;
display y1_4;
display z1_4;
display x1_5;
display y1_5;
display z1_5;
display x1_6;
display y1_6;
display z1_6;
display x1_7;
display y1_7;
display z1_7;
display x1_8;
display y1_8;
display z1_8;
display x1_9;
display y1_9;
display z1_9;
display x1_10;
display y1_10;
display z1_10;
display x2_1;
display y2_1;
display z2_1;
display x2_2;
display y2_2;
display z2_2;
display x2_3;
display y2_3;
display z2_3;
display x2_4;
display y2_4;
display z2_4;
display x2_5;
display y2_5;
display z2_5;
display x2_6;
display y2_6;
display z2_6;
display x2_7;
display y2_7;
display z2_7;
display x2_8;
display y2_8;
display z2_8;
display x2_9;
display y2_9;
display z2_9;
display x2_10;
display y2_10;
display z2_10;
display x3_1;
display y3_1;
display z3_1;
display x3_2;
display y3_2;
display z3_2;
display x3_3;
display y3_3;
display z3_3;
display x3_4;
display y3_4;
display z3_4;
display x3_5;
display y3_5;
display z3_5;
display x3_6;
display y3_6;
display z3_6;
display x3_7;
display y3_7;
display z3_7;
display x3_8;
display y3_8;
display z3_8;
display x3_9;
display y3_9;
display z3_9;
display x3_10;
display y3_10;
display z3_10;
display x4_1;
display y4_1;
display z4_1;
display x4_2;
display y4_2;
display z4_2;
display x4_3;
display y4_3;
display z4_3;
display x4_4;
display y4_4;
display z4_4;
display x4_5;
display y4_5;
display z4_5;
display x4_6;
display y4_6;
display z4_6;
display x4_7;
display y4_7;
display z4_7;
display x4_8;
display y4_8;
display z4_8;
display x4_9;
display y4_9;
display z4_9;
display x4_10;
display y4_10;
display z4_10;
display x5_1;
display y5_1;
display z5_1;
display x5_2;
display y5_2;
display z5_2;
display x5_3;
display y5_3;
display z5_3;
display x5_4;
display y5_4;
display z5_4;
display x5_5;
display y5_5;
display z5_5;
display x5_6;
display y5_6;
display z5_6;
display x5_7;
display y5_7;
display z5_7;
display x5_8;
display y5_8;
display z5_8;
display x5_9;
display y5_9;
display z5_9;
display x5_10;
display y5_10;
display z5_10;
display x6_1;
display y6_1;
display z6_1;
display x6_2;
display y6_2;
display z6_2;
display x6_3;
display y6_3;
display z6_3;
display x6_4;
display y6_4;
display z6_4;
display x6_5;
display y6_5;
display z6_5;
display x6_6;
display y6_6;
display z6_6;
display x6_7;
display y6_7;
display z6_7;
display x6_8;
display y6_8;
display z6_8;
display x6_9;
display y6_9;
display z6_9;
display x6_10;
display y6_10;
display z6_10;
display x7_1;
display y7_1;
display z7_1;
display x7_2;
display y7_2;
display z7_2;
display x7_3;
display y7_3;
display z7_3;
display x7_4;
display y7_4;
display z7_4;
display x7_5;
display y7_5;
display z7_5;
display x7_6;
display y7_6;
display z7_6;
display x7_7;
display y7_7;
display z7_7;
display x7_8;
display y7_8;
display z7_8;
display x7_9;
display y7_9;
display z7_9;
display x7_10;
display y7_10;
display z7_10;
display x8_1;
display y8_1;
display z8_1;
display x8_2;
display y8_2;
display z8_2;
display x8_3;
display y8_3;
display z8_3;
display x8_4;
display y8_4;
display z8_4;
display x8_5;
display y8_5;
display z8_5;
display x8_6;
display y8_6;
display z8_6;
display x8_7;
display y8_7;
display z8_7;
display x8_8;
display y8_8;
display z8_8;
display x8_9;
display y8_9;
display z8_9;
display x8_10;
display y8_10;
display z8_10;
display x9_1;
display y9_1;
display z9_1;
display x9_2;
display y9_2;
display z9_2;
display x9_3;
display y9_3;
display z9_3;
display x9_4;
display y9_4;
display z9_4;
display x9_5;
display y9_5;
display z9_5;
display x9_6;
display y9_6;
display z9_6;
display x9_7;
display y9_7;
display z9_7;
display x9_8;
display y9_8;
display z9_8;
display x9_9;
display y9_9;
display z9_9;
display x9_10;
display y9_10;
display z9_10;
display x10_1;
display y10_1;
display z10_1;
display x10_2;
display y10_2;
display z10_2;
display x10_3;
display y10_3;
display z10_3;
display x10_4;
display y10_4;
display z10_4;
display x10_5;
display y10_5;
display z10_5;
display x10_6;
display y10_6;
display z10_6;
display x10_7;
display y10_7;
display z10_7;
display x10_8;
display y10_8;
display z10_8;
display x10_9;
display y10_9;
display z10_9;
display x10_10;
display y10_10;
display z10_10;
display obj;
