#
#**************************
# SET UP THE INITIAL DATA *
#**************************
#   Problem :
#   *********
#   Find a polygon of maximal area, under the constraint that no two of
#   its vertices are further apart than 1.  The considered class of polygons
#   has one vertex fixed at the origin and one of its neighbouring vertices
#   constrained to lie on the horizontal axis.  Furthermore, vertices
#   must have increasing abscissas, which implies that the polygons are
#   of the form
#
#                     ---*-----*-*
#                    *            \
#                   /              *
#                  /                \
#                  *-----------------*
#                (0,0)
#   Source:
#   Ph. Toint
#   SIF input: Ph. Toint, June 1990.
#   classification QQR2-AY-V-V
#   HNS is the number of sides - 1
#   The number of variables in the problem is N = 2 * ( HNS + 1 )
#   of which 3 are fixed and one bounded above.
#IE HNS                 3              $ n = 8
#IE HNS                 5              $ n = 12
#IE HNS                 24             $ n = 50
#IE HNS                 249            $ n = 500
#   Constants
#   Computed parameters
#   Objective function
#   Ordering constraints on the abscissas
#   Constraints on the distance between two vertices
#   Objective function
#   Constraints on the distances between vertices
#   Objective function
#   Constraints on the distances between intermediate points
#   Lower bound given by the area of the intersection of the unit sphere
#   and the positive orthant
#   Solution
#LO SOLTN(3)            -0.5
#LO SOLTN(5)            -0.6203656662
#LO SOLTN(25)           -0.6410147561
	param hns := 49;
	param hnsm1 := -1 + (49);
	param ip1 := 1 + (48);
	param im1 := -1 + (49);

	var x0 >= 0.0 ,  <= 0.0;
	var y0 >= 0.0 ,  <= 0.0;
	var x1 >= 0.0;
	var y1 >= 0.0;
	var x2 >= 0.0;
	var y2 >= 0.0;
	var x3 >= 0.0;
	var y3 >= 0.0;
	var x4 >= 0.0;
	var y4 >= 0.0;
	var x5 >= 0.0;
	var y5 >= 0.0;
	var x6 >= 0.0;
	var y6 >= 0.0;
	var x7 >= 0.0;
	var y7 >= 0.0;
	var x8 >= 0.0;
	var y8 >= 0.0;
	var x9 >= 0.0;
	var y9 >= 0.0;
	var x10 >= 0.0;
	var y10 >= 0.0;
	var x11 >= 0.0;
	var y11 >= 0.0;
	var x12 >= 0.0;
	var y12 >= 0.0;
	var x13 >= 0.0;
	var y13 >= 0.0;
	var x14 >= 0.0;
	var y14 >= 0.0;
	var x15 >= 0.0;
	var y15 >= 0.0;
	var x16 >= 0.0;
	var y16 >= 0.0;
	var x17 >= 0.0;
	var y17 >= 0.0;
	var x18 >= 0.0;
	var y18 >= 0.0;
	var x19 >= 0.0;
	var y19 >= 0.0;
	var x20 >= 0.0;
	var y20 >= 0.0;
	var x21 >= 0.0;
	var y21 >= 0.0;
	var x22 >= 0.0;
	var y22 >= 0.0;
	var x23 >= 0.0;
	var y23 >= 0.0;
	var x24 >= 0.0;
	var y24 >= 0.0;
	var x25 >= 0.0;
	var y25 >= 0.0;
	var x26 >= 0.0;
	var y26 >= 0.0;
	var x27 >= 0.0;
	var y27 >= 0.0;
	var x28 >= 0.0;
	var y28 >= 0.0;
	var x29 >= 0.0;
	var y29 >= 0.0;
	var x30 >= 0.0;
	var y30 >= 0.0;
	var x31 >= 0.0;
	var y31 >= 0.0;
	var x32 >= 0.0;
	var y32 >= 0.0;
	var x33 >= 0.0;
	var y33 >= 0.0;
	var x34 >= 0.0;
	var y34 >= 0.0;
	var x35 >= 0.0;
	var y35 >= 0.0;
	var x36 >= 0.0;
	var y36 >= 0.0;
	var x37 >= 0.0;
	var y37 >= 0.0;
	var x38 >= 0.0;
	var y38 >= 0.0;
	var x39 >= 0.0;
	var y39 >= 0.0;
	var x40 >= 0.0;
	var y40 >= 0.0;
	var x41 >= 0.0;
	var y41 >= 0.0;
	var x42 >= 0.0;
	var y42 >= 0.0;
	var x43 >= 0.0;
	var y43 >= 0.0;
	var x44 >= 0.0;
	var y44 >= 0.0;
	var x45 >= 0.0;
	var y45 >= 0.0;
	var x46 >= 0.0;
	var y46 >= 0.0;
	var x47 >= 0.0;
	var y47 >= 0.0;
	var x48 >= 0.0;
	var y48 >= 0.0;
	var x49 >= 0.0 ,  <= 1.0 ,  := 2.0;
	var y49 >= 0.0 ,  <= 0.0;

minimize obj:
	 - 0.5*(x1-x0) * (y1+y0) - 0.5*(x2-x1) * (y2+y1) - 0.5*(x3-x2) * (y3+y2) - 
	0.5*(x4-x3) * (y4+y3) - 0.5*(x5-x4) * (y5+y4) - 0.5*(x6-x5) * (y6+y5) - 
	0.5*(x7-x6) * (y7+y6) - 0.5*(x8-x7) * (y8+y7) - 0.5*(x9-x8) * (y9+y8) - 
	0.5*(x10-x9) * (y10+y9) - 0.5*(x11-x10) * (y11+y10) - 0.5*(x12-x11) * (y12+y11) 
	- 0.5*(x13-x12) * (y13+y12) - 0.5*(x14-x13) * (y14+y13) - 0.5*(x15-x14) * 
	(y15+y14) - 0.5*(x16-x15) * (y16+y15) - 0.5*(x17-x16) * (y17+y16) - 
	0.5*(x18-x17) * (y18+y17) - 0.5*(x19-x18) * (y19+y18) - 0.5*(x20-x19) * 
	(y20+y19) - 0.5*(x21-x20) * (y21+y20) - 0.5*(x22-x21) * (y22+y21) - 
	0.5*(x23-x22) * (y23+y22) - 0.5*(x24-x23) * (y24+y23) - 0.5*(x25-x24) * 
	(y25+y24) - 0.5*(x26-x25) * (y26+y25) - 0.5*(x27-x26) * (y27+y26) - 
	0.5*(x28-x27) * (y28+y27) - 0.5*(x29-x28) * (y29+y28) - 0.5*(x30-x29) * 
	(y30+y29) - 0.5*(x31-x30) * (y31+y30) - 0.5*(x32-x31) * (y32+y31) - 
	0.5*(x33-x32) * (y33+y32) - 0.5*(x34-x33) * (y34+y33) - 0.5*(x35-x34) * 
	(y35+y34) - 0.5*(x36-x35) * (y36+y35) - 0.5*(x37-x36) * (y37+y36) - 
	0.5*(x38-x37) * (y38+y37) - 0.5*(x39-x38) * (y39+y38) - 0.5*(x40-x39) * 
	(y40+y39) - 0.5*(x41-x40) * (y41+y40) - 0.5*(x42-x41) * (y42+y41) - 
	0.5*(x43-x42) * (y43+y42) - 0.5*(x44-x43) * (y44+y43) - 0.5*(x45-x44) * 
	(y45+y44) - 0.5*(x46-x45) * (y46+y45) - 0.5*(x47-x46) * (y47+y46) - 
	0.5*(x48-x47) * (y48+y47) - 0.5*(x49-x48) * (y49+y48);

subject to ox1:
	0 <= x2 - x1;
subject to ox2:
	0 <= x3 - x2;
subject to ox3:
	0 <= x4 - x3;
subject to ox4:
	0 <= x5 - x4;
subject to ox5:
	0 <= x6 - x5;
subject to ox6:
	0 <= x7 - x6;
subject to ox7:
	0 <= x8 - x7;
subject to ox8:
	0 <= x9 - x8;
subject to ox9:
	0 <= x10 - x9;
subject to ox10:
	0 <= x11 - x10;
subject to ox11:
	0 <= x12 - x11;
subject to ox12:
	0 <= x13 - x12;
subject to ox13:
	0 <= x14 - x13;
subject to ox14:
	0 <= x15 - x14;
subject to ox15:
	0 <= x16 - x15;
subject to ox16:
	0 <= x17 - x16;
subject to ox17:
	0 <= x18 - x17;
subject to ox18:
	0 <= x19 - x18;
subject to ox19:
	0 <= x20 - x19;
subject to ox20:
	0 <= x21 - x20;
subject to ox21:
	0 <= x22 - x21;
subject to ox22:
	0 <= x23 - x22;
subject to ox23:
	0 <= x24 - x23;
subject to ox24:
	0 <= x25 - x24;
subject to ox25:
	0 <= x26 - x25;
subject to ox26:
	0 <= x27 - x26;
subject to ox27:
	0 <= x28 - x27;
subject to ox28:
	0 <= x29 - x28;
subject to ox29:
	0 <= x30 - x29;
subject to ox30:
	0 <= x31 - x30;
subject to ox31:
	0 <= x32 - x31;
subject to ox32:
	0 <= x33 - x32;
subject to ox33:
	0 <= x34 - x33;
subject to ox34:
	0 <= x35 - x34;
subject to ox35:
	0 <= x36 - x35;
subject to ox36:
	0 <= x37 - x36;
subject to ox37:
	0 <= x38 - x37;
subject to ox38:
	0 <= x39 - x38;
subject to ox39:
	0 <= x40 - x39;
subject to ox40:
	0 <= x41 - x40;
subject to ox41:
	0 <= x42 - x41;
subject to ox42:
	0 <= x43 - x42;
subject to ox43:
	0 <= x44 - x43;
subject to ox44:
	0 <= x45 - x44;
subject to ox45:
	0 <= x46 - x45;
subject to ox46:
	0 <= x47 - x46;
subject to ox47:
	0 <= x48 - x47;
subject to ox48:
	0 <= x49 - x48;
subject to d1_0:
	0 >= (x1-x0) * (x1-x0) + (y1-y0) * (y1-y0) - 1.0;
subject to d2_0:
	0 >= (x2-x0) * (x2-x0) + (y2-y0) * (y2-y0) - 1.0;
subject to d2_1:
	0 >= (x2-x1) * (x2-x1) + (y2-y1) * (y2-y1) - 1.0;
subject to d3_0:
	0 >= (x3-x0) * (x3-x0) + (y3-y0) * (y3-y0) - 1.0;
subject to d3_1:
	0 >= (x3-x1) * (x3-x1) + (y3-y1) * (y3-y1) - 1.0;
subject to d3_2:
	0 >= (x3-x2) * (x3-x2) + (y3-y2) * (y3-y2) - 1.0;
subject to d4_0:
	0 >= (x4-x0) * (x4-x0) + (y4-y0) * (y4-y0) - 1.0;
subject to d4_1:
	0 >= (x4-x1) * (x4-x1) + (y4-y1) * (y4-y1) - 1.0;
subject to d4_2:
	0 >= (x4-x2) * (x4-x2) + (y4-y2) * (y4-y2) - 1.0;
subject to d4_3:
	0 >= (x4-x3) * (x4-x3) + (y4-y3) * (y4-y3) - 1.0;
subject to d5_0:
	0 >= (x5-x0) * (x5-x0) + (y5-y0) * (y5-y0) - 1.0;
subject to d5_1:
	0 >= (x5-x1) * (x5-x1) + (y5-y1) * (y5-y1) - 1.0;
subject to d5_2:
	0 >= (x5-x2) * (x5-x2) + (y5-y2) * (y5-y2) - 1.0;
subject to d5_3:
	0 >= (x5-x3) * (x5-x3) + (y5-y3) * (y5-y3) - 1.0;
subject to d5_4:
	0 >= (x5-x4) * (x5-x4) + (y5-y4) * (y5-y4) - 1.0;
subject to d6_0:
	0 >= (x6-x0) * (x6-x0) + (y6-y0) * (y6-y0) - 1.0;
subject to d6_1:
	0 >= (x6-x1) * (x6-x1) + (y6-y1) * (y6-y1) - 1.0;
subject to d6_2:
	0 >= (x6-x2) * (x6-x2) + (y6-y2) * (y6-y2) - 1.0;
subject to d6_3:
	0 >= (x6-x3) * (x6-x3) + (y6-y3) * (y6-y3) - 1.0;
subject to d6_4:
	0 >= (x6-x4) * (x6-x4) + (y6-y4) * (y6-y4) - 1.0;
subject to d6_5:
	0 >= (x6-x5) * (x6-x5) + (y6-y5) * (y6-y5) - 1.0;
subject to d7_0:
	0 >= (x7-x0) * (x7-x0) + (y7-y0) * (y7-y0) - 1.0;
subject to d7_1:
	0 >= (x7-x1) * (x7-x1) + (y7-y1) * (y7-y1) - 1.0;
subject to d7_2:
	0 >= (x7-x2) * (x7-x2) + (y7-y2) * (y7-y2) - 1.0;
subject to d7_3:
	0 >= (x7-x3) * (x7-x3) + (y7-y3) * (y7-y3) - 1.0;
subject to d7_4:
	0 >= (x7-x4) * (x7-x4) + (y7-y4) * (y7-y4) - 1.0;
subject to d7_5:
	0 >= (x7-x5) * (x7-x5) + (y7-y5) * (y7-y5) - 1.0;
subject to d7_6:
	0 >= (x7-x6) * (x7-x6) + (y7-y6) * (y7-y6) - 1.0;
subject to d8_0:
	0 >= (x8-x0) * (x8-x0) + (y8-y0) * (y8-y0) - 1.0;
subject to d8_1:
	0 >= (x8-x1) * (x8-x1) + (y8-y1) * (y8-y1) - 1.0;
subject to d8_2:
	0 >= (x8-x2) * (x8-x2) + (y8-y2) * (y8-y2) - 1.0;
subject to d8_3:
	0 >= (x8-x3) * (x8-x3) + (y8-y3) * (y8-y3) - 1.0;
subject to d8_4:
	0 >= (x8-x4) * (x8-x4) + (y8-y4) * (y8-y4) - 1.0;
subject to d8_5:
	0 >= (x8-x5) * (x8-x5) + (y8-y5) * (y8-y5) - 1.0;
subject to d8_6:
	0 >= (x8-x6) * (x8-x6) + (y8-y6) * (y8-y6) - 1.0;
subject to d8_7:
	0 >= (x8-x7) * (x8-x7) + (y8-y7) * (y8-y7) - 1.0;
subject to d9_0:
	0 >= (x9-x0) * (x9-x0) + (y9-y0) * (y9-y0) - 1.0;
subject to d9_1:
	0 >= (x9-x1) * (x9-x1) + (y9-y1) * (y9-y1) - 1.0;
subject to d9_2:
	0 >= (x9-x2) * (x9-x2) + (y9-y2) * (y9-y2) - 1.0;
subject to d9_3:
	0 >= (x9-x3) * (x9-x3) + (y9-y3) * (y9-y3) - 1.0;
subject to d9_4:
	0 >= (x9-x4) * (x9-x4) + (y9-y4) * (y9-y4) - 1.0;
subject to d9_5:
	0 >= (x9-x5) * (x9-x5) + (y9-y5) * (y9-y5) - 1.0;
subject to d9_6:
	0 >= (x9-x6) * (x9-x6) + (y9-y6) * (y9-y6) - 1.0;
subject to d9_7:
	0 >= (x9-x7) * (x9-x7) + (y9-y7) * (y9-y7) - 1.0;
subject to d9_8:
	0 >= (x9-x8) * (x9-x8) + (y9-y8) * (y9-y8) - 1.0;
subject to d10_0:
	0 >= (x10-x0) * (x10-x0) + (y10-y0) * (y10-y0) - 1.0;
subject to d10_1:
	0 >= (x10-x1) * (x10-x1) + (y10-y1) * (y10-y1) - 1.0;
subject to d10_2:
	0 >= (x10-x2) * (x10-x2) + (y10-y2) * (y10-y2) - 1.0;
subject to d10_3:
	0 >= (x10-x3) * (x10-x3) + (y10-y3) * (y10-y3) - 1.0;
subject to d10_4:
	0 >= (x10-x4) * (x10-x4) + (y10-y4) * (y10-y4) - 1.0;
subject to d10_5:
	0 >= (x10-x5) * (x10-x5) + (y10-y5) * (y10-y5) - 1.0;
subject to d10_6:
	0 >= (x10-x6) * (x10-x6) + (y10-y6) * (y10-y6) - 1.0;
subject to d10_7:
	0 >= (x10-x7) * (x10-x7) + (y10-y7) * (y10-y7) - 1.0;
subject to d10_8:
	0 >= (x10-x8) * (x10-x8) + (y10-y8) * (y10-y8) - 1.0;
subject to d10_9:
	0 >= (x10-x9) * (x10-x9) + (y10-y9) * (y10-y9) - 1.0;
subject to d11_0:
	0 >= (x11-x0) * (x11-x0) + (y11-y0) * (y11-y0) - 1.0;
subject to d11_1:
	0 >= (x11-x1) * (x11-x1) + (y11-y1) * (y11-y1) - 1.0;
subject to d11_2:
	0 >= (x11-x2) * (x11-x2) + (y11-y2) * (y11-y2) - 1.0;
subject to d11_3:
	0 >= (x11-x3) * (x11-x3) + (y11-y3) * (y11-y3) - 1.0;
subject to d11_4:
	0 >= (x11-x4) * (x11-x4) + (y11-y4) * (y11-y4) - 1.0;
subject to d11_5:
	0 >= (x11-x5) * (x11-x5) + (y11-y5) * (y11-y5) - 1.0;
subject to d11_6:
	0 >= (x11-x6) * (x11-x6) + (y11-y6) * (y11-y6) - 1.0;
subject to d11_7:
	0 >= (x11-x7) * (x11-x7) + (y11-y7) * (y11-y7) - 1.0;
subject to d11_8:
	0 >= (x11-x8) * (x11-x8) + (y11-y8) * (y11-y8) - 1.0;
subject to d11_9:
	0 >= (x11-x9) * (x11-x9) + (y11-y9) * (y11-y9) - 1.0;
subject to d11_10:
	0 >= (x11-x10) * (x11-x10) + (y11-y10) * (y11-y10) - 1.0;
subject to d12_0:
	0 >= (x12-x0) * (x12-x0) + (y12-y0) * (y12-y0) - 1.0;
subject to d12_1:
	0 >= (x12-x1) * (x12-x1) + (y12-y1) * (y12-y1) - 1.0;
subject to d12_2:
	0 >= (x12-x2) * (x12-x2) + (y12-y2) * (y12-y2) - 1.0;
subject to d12_3:
	0 >= (x12-x3) * (x12-x3) + (y12-y3) * (y12-y3) - 1.0;
subject to d12_4:
	0 >= (x12-x4) * (x12-x4) + (y12-y4) * (y12-y4) - 1.0;
subject to d12_5:
	0 >= (x12-x5) * (x12-x5) + (y12-y5) * (y12-y5) - 1.0;
subject to d12_6:
	0 >= (x12-x6) * (x12-x6) + (y12-y6) * (y12-y6) - 1.0;
subject to d12_7:
	0 >= (x12-x7) * (x12-x7) + (y12-y7) * (y12-y7) - 1.0;
subject to d12_8:
	0 >= (x12-x8) * (x12-x8) + (y12-y8) * (y12-y8) - 1.0;
subject to d12_9:
	0 >= (x12-x9) * (x12-x9) + (y12-y9) * (y12-y9) - 1.0;
subject to d12_10:
	0 >= (x12-x10) * (x12-x10) + (y12-y10) * (y12-y10) - 1.0;
subject to d12_11:
	0 >= (x12-x11) * (x12-x11) + (y12-y11) * (y12-y11) - 1.0;
subject to d13_0:
	0 >= (x13-x0) * (x13-x0) + (y13-y0) * (y13-y0) - 1.0;
subject to d13_1:
	0 >= (x13-x1) * (x13-x1) + (y13-y1) * (y13-y1) - 1.0;
subject to d13_2:
	0 >= (x13-x2) * (x13-x2) + (y13-y2) * (y13-y2) - 1.0;
subject to d13_3:
	0 >= (x13-x3) * (x13-x3) + (y13-y3) * (y13-y3) - 1.0;
subject to d13_4:
	0 >= (x13-x4) * (x13-x4) + (y13-y4) * (y13-y4) - 1.0;
subject to d13_5:
	0 >= (x13-x5) * (x13-x5) + (y13-y5) * (y13-y5) - 1.0;
subject to d13_6:
	0 >= (x13-x6) * (x13-x6) + (y13-y6) * (y13-y6) - 1.0;
subject to d13_7:
	0 >= (x13-x7) * (x13-x7) + (y13-y7) * (y13-y7) - 1.0;
subject to d13_8:
	0 >= (x13-x8) * (x13-x8) + (y13-y8) * (y13-y8) - 1.0;
subject to d13_9:
	0 >= (x13-x9) * (x13-x9) + (y13-y9) * (y13-y9) - 1.0;
subject to d13_10:
	0 >= (x13-x10) * (x13-x10) + (y13-y10) * (y13-y10) - 1.0;
subject to d13_11:
	0 >= (x13-x11) * (x13-x11) + (y13-y11) * (y13-y11) - 1.0;
subject to d13_12:
	0 >= (x13-x12) * (x13-x12) + (y13-y12) * (y13-y12) - 1.0;
subject to d14_0:
	0 >= (x14-x0) * (x14-x0) + (y14-y0) * (y14-y0) - 1.0;
subject to d14_1:
	0 >= (x14-x1) * (x14-x1) + (y14-y1) * (y14-y1) - 1.0;
subject to d14_2:
	0 >= (x14-x2) * (x14-x2) + (y14-y2) * (y14-y2) - 1.0;
subject to d14_3:
	0 >= (x14-x3) * (x14-x3) + (y14-y3) * (y14-y3) - 1.0;
subject to d14_4:
	0 >= (x14-x4) * (x14-x4) + (y14-y4) * (y14-y4) - 1.0;
subject to d14_5:
	0 >= (x14-x5) * (x14-x5) + (y14-y5) * (y14-y5) - 1.0;
subject to d14_6:
	0 >= (x14-x6) * (x14-x6) + (y14-y6) * (y14-y6) - 1.0;
subject to d14_7:
	0 >= (x14-x7) * (x14-x7) + (y14-y7) * (y14-y7) - 1.0;
subject to d14_8:
	0 >= (x14-x8) * (x14-x8) + (y14-y8) * (y14-y8) - 1.0;
subject to d14_9:
	0 >= (x14-x9) * (x14-x9) + (y14-y9) * (y14-y9) - 1.0;
subject to d14_10:
	0 >= (x14-x10) * (x14-x10) + (y14-y10) * (y14-y10) - 1.0;
subject to d14_11:
	0 >= (x14-x11) * (x14-x11) + (y14-y11) * (y14-y11) - 1.0;
subject to d14_12:
	0 >= (x14-x12) * (x14-x12) + (y14-y12) * (y14-y12) - 1.0;
subject to d14_13:
	0 >= (x14-x13) * (x14-x13) + (y14-y13) * (y14-y13) - 1.0;
subject to d15_0:
	0 >= (x15-x0) * (x15-x0) + (y15-y0) * (y15-y0) - 1.0;
subject to d15_1:
	0 >= (x15-x1) * (x15-x1) + (y15-y1) * (y15-y1) - 1.0;
subject to d15_2:
	0 >= (x15-x2) * (x15-x2) + (y15-y2) * (y15-y2) - 1.0;
subject to d15_3:
	0 >= (x15-x3) * (x15-x3) + (y15-y3) * (y15-y3) - 1.0;
subject to d15_4:
	0 >= (x15-x4) * (x15-x4) + (y15-y4) * (y15-y4) - 1.0;
subject to d15_5:
	0 >= (x15-x5) * (x15-x5) + (y15-y5) * (y15-y5) - 1.0;
subject to d15_6:
	0 >= (x15-x6) * (x15-x6) + (y15-y6) * (y15-y6) - 1.0;
subject to d15_7:
	0 >= (x15-x7) * (x15-x7) + (y15-y7) * (y15-y7) - 1.0;
subject to d15_8:
	0 >= (x15-x8) * (x15-x8) + (y15-y8) * (y15-y8) - 1.0;
subject to d15_9:
	0 >= (x15-x9) * (x15-x9) + (y15-y9) * (y15-y9) - 1.0;
subject to d15_10:
	0 >= (x15-x10) * (x15-x10) + (y15-y10) * (y15-y10) - 1.0;
subject to d15_11:
	0 >= (x15-x11) * (x15-x11) + (y15-y11) * (y15-y11) - 1.0;
subject to d15_12:
	0 >= (x15-x12) * (x15-x12) + (y15-y12) * (y15-y12) - 1.0;
subject to d15_13:
	0 >= (x15-x13) * (x15-x13) + (y15-y13) * (y15-y13) - 1.0;
subject to d15_14:
	0 >= (x15-x14) * (x15-x14) + (y15-y14) * (y15-y14) - 1.0;
subject to d16_0:
	0 >= (x16-x0) * (x16-x0) + (y16-y0) * (y16-y0) - 1.0;
subject to d16_1:
	0 >= (x16-x1) * (x16-x1) + (y16-y1) * (y16-y1) - 1.0;
subject to d16_2:
	0 >= (x16-x2) * (x16-x2) + (y16-y2) * (y16-y2) - 1.0;
subject to d16_3:
	0 >= (x16-x3) * (x16-x3) + (y16-y3) * (y16-y3) - 1.0;
subject to d16_4:
	0 >= (x16-x4) * (x16-x4) + (y16-y4) * (y16-y4) - 1.0;
subject to d16_5:
	0 >= (x16-x5) * (x16-x5) + (y16-y5) * (y16-y5) - 1.0;
subject to d16_6:
	0 >= (x16-x6) * (x16-x6) + (y16-y6) * (y16-y6) - 1.0;
subject to d16_7:
	0 >= (x16-x7) * (x16-x7) + (y16-y7) * (y16-y7) - 1.0;
subject to d16_8:
	0 >= (x16-x8) * (x16-x8) + (y16-y8) * (y16-y8) - 1.0;
subject to d16_9:
	0 >= (x16-x9) * (x16-x9) + (y16-y9) * (y16-y9) - 1.0;
subject to d16_10:
	0 >= (x16-x10) * (x16-x10) + (y16-y10) * (y16-y10) - 1.0;
subject to d16_11:
	0 >= (x16-x11) * (x16-x11) + (y16-y11) * (y16-y11) - 1.0;
subject to d16_12:
	0 >= (x16-x12) * (x16-x12) + (y16-y12) * (y16-y12) - 1.0;
subject to d16_13:
	0 >= (x16-x13) * (x16-x13) + (y16-y13) * (y16-y13) - 1.0;
subject to d16_14:
	0 >= (x16-x14) * (x16-x14) + (y16-y14) * (y16-y14) - 1.0;
subject to d16_15:
	0 >= (x16-x15) * (x16-x15) + (y16-y15) * (y16-y15) - 1.0;
subject to d17_0:
	0 >= (x17-x0) * (x17-x0) + (y17-y0) * (y17-y0) - 1.0;
subject to d17_1:
	0 >= (x17-x1) * (x17-x1) + (y17-y1) * (y17-y1) - 1.0;
subject to d17_2:
	0 >= (x17-x2) * (x17-x2) + (y17-y2) * (y17-y2) - 1.0;
subject to d17_3:
	0 >= (x17-x3) * (x17-x3) + (y17-y3) * (y17-y3) - 1.0;
subject to d17_4:
	0 >= (x17-x4) * (x17-x4) + (y17-y4) * (y17-y4) - 1.0;
subject to d17_5:
	0 >= (x17-x5) * (x17-x5) + (y17-y5) * (y17-y5) - 1.0;
subject to d17_6:
	0 >= (x17-x6) * (x17-x6) + (y17-y6) * (y17-y6) - 1.0;
subject to d17_7:
	0 >= (x17-x7) * (x17-x7) + (y17-y7) * (y17-y7) - 1.0;
subject to d17_8:
	0 >= (x17-x8) * (x17-x8) + (y17-y8) * (y17-y8) - 1.0;
subject to d17_9:
	0 >= (x17-x9) * (x17-x9) + (y17-y9) * (y17-y9) - 1.0;
subject to d17_10:
	0 >= (x17-x10) * (x17-x10) + (y17-y10) * (y17-y10) - 1.0;
subject to d17_11:
	0 >= (x17-x11) * (x17-x11) + (y17-y11) * (y17-y11) - 1.0;
subject to d17_12:
	0 >= (x17-x12) * (x17-x12) + (y17-y12) * (y17-y12) - 1.0;
subject to d17_13:
	0 >= (x17-x13) * (x17-x13) + (y17-y13) * (y17-y13) - 1.0;
subject to d17_14:
	0 >= (x17-x14) * (x17-x14) + (y17-y14) * (y17-y14) - 1.0;
subject to d17_15:
	0 >= (x17-x15) * (x17-x15) + (y17-y15) * (y17-y15) - 1.0;
subject to d17_16:
	0 >= (x17-x16) * (x17-x16) + (y17-y16) * (y17-y16) - 1.0;
subject to d18_0:
	0 >= (x18-x0) * (x18-x0) + (y18-y0) * (y18-y0) - 1.0;
subject to d18_1:
	0 >= (x18-x1) * (x18-x1) + (y18-y1) * (y18-y1) - 1.0;
subject to d18_2:
	0 >= (x18-x2) * (x18-x2) + (y18-y2) * (y18-y2) - 1.0;
subject to d18_3:
	0 >= (x18-x3) * (x18-x3) + (y18-y3) * (y18-y3) - 1.0;
subject to d18_4:
	0 >= (x18-x4) * (x18-x4) + (y18-y4) * (y18-y4) - 1.0;
subject to d18_5:
	0 >= (x18-x5) * (x18-x5) + (y18-y5) * (y18-y5) - 1.0;
subject to d18_6:
	0 >= (x18-x6) * (x18-x6) + (y18-y6) * (y18-y6) - 1.0;
subject to d18_7:
	0 >= (x18-x7) * (x18-x7) + (y18-y7) * (y18-y7) - 1.0;
subject to d18_8:
	0 >= (x18-x8) * (x18-x8) + (y18-y8) * (y18-y8) - 1.0;
subject to d18_9:
	0 >= (x18-x9) * (x18-x9) + (y18-y9) * (y18-y9) - 1.0;
subject to d18_10:
	0 >= (x18-x10) * (x18-x10) + (y18-y10) * (y18-y10) - 1.0;
subject to d18_11:
	0 >= (x18-x11) * (x18-x11) + (y18-y11) * (y18-y11) - 1.0;
subject to d18_12:
	0 >= (x18-x12) * (x18-x12) + (y18-y12) * (y18-y12) - 1.0;
subject to d18_13:
	0 >= (x18-x13) * (x18-x13) + (y18-y13) * (y18-y13) - 1.0;
subject to d18_14:
	0 >= (x18-x14) * (x18-x14) + (y18-y14) * (y18-y14) - 1.0;
subject to d18_15:
	0 >= (x18-x15) * (x18-x15) + (y18-y15) * (y18-y15) - 1.0;
subject to d18_16:
	0 >= (x18-x16) * (x18-x16) + (y18-y16) * (y18-y16) - 1.0;
subject to d18_17:
	0 >= (x18-x17) * (x18-x17) + (y18-y17) * (y18-y17) - 1.0;
subject to d19_0:
	0 >= (x19-x0) * (x19-x0) + (y19-y0) * (y19-y0) - 1.0;
subject to d19_1:
	0 >= (x19-x1) * (x19-x1) + (y19-y1) * (y19-y1) - 1.0;
subject to d19_2:
	0 >= (x19-x2) * (x19-x2) + (y19-y2) * (y19-y2) - 1.0;
subject to d19_3:
	0 >= (x19-x3) * (x19-x3) + (y19-y3) * (y19-y3) - 1.0;
subject to d19_4:
	0 >= (x19-x4) * (x19-x4) + (y19-y4) * (y19-y4) - 1.0;
subject to d19_5:
	0 >= (x19-x5) * (x19-x5) + (y19-y5) * (y19-y5) - 1.0;
subject to d19_6:
	0 >= (x19-x6) * (x19-x6) + (y19-y6) * (y19-y6) - 1.0;
subject to d19_7:
	0 >= (x19-x7) * (x19-x7) + (y19-y7) * (y19-y7) - 1.0;
subject to d19_8:
	0 >= (x19-x8) * (x19-x8) + (y19-y8) * (y19-y8) - 1.0;
subject to d19_9:
	0 >= (x19-x9) * (x19-x9) + (y19-y9) * (y19-y9) - 1.0;
subject to d19_10:
	0 >= (x19-x10) * (x19-x10) + (y19-y10) * (y19-y10) - 1.0;
subject to d19_11:
	0 >= (x19-x11) * (x19-x11) + (y19-y11) * (y19-y11) - 1.0;
subject to d19_12:
	0 >= (x19-x12) * (x19-x12) + (y19-y12) * (y19-y12) - 1.0;
subject to d19_13:
	0 >= (x19-x13) * (x19-x13) + (y19-y13) * (y19-y13) - 1.0;
subject to d19_14:
	0 >= (x19-x14) * (x19-x14) + (y19-y14) * (y19-y14) - 1.0;
subject to d19_15:
	0 >= (x19-x15) * (x19-x15) + (y19-y15) * (y19-y15) - 1.0;
subject to d19_16:
	0 >= (x19-x16) * (x19-x16) + (y19-y16) * (y19-y16) - 1.0;
subject to d19_17:
	0 >= (x19-x17) * (x19-x17) + (y19-y17) * (y19-y17) - 1.0;
subject to d19_18:
	0 >= (x19-x18) * (x19-x18) + (y19-y18) * (y19-y18) - 1.0;
subject to d20_0:
	0 >= (x20-x0) * (x20-x0) + (y20-y0) * (y20-y0) - 1.0;
subject to d20_1:
	0 >= (x20-x1) * (x20-x1) + (y20-y1) * (y20-y1) - 1.0;
subject to d20_2:
	0 >= (x20-x2) * (x20-x2) + (y20-y2) * (y20-y2) - 1.0;
subject to d20_3:
	0 >= (x20-x3) * (x20-x3) + (y20-y3) * (y20-y3) - 1.0;
subject to d20_4:
	0 >= (x20-x4) * (x20-x4) + (y20-y4) * (y20-y4) - 1.0;
subject to d20_5:
	0 >= (x20-x5) * (x20-x5) + (y20-y5) * (y20-y5) - 1.0;
subject to d20_6:
	0 >= (x20-x6) * (x20-x6) + (y20-y6) * (y20-y6) - 1.0;
subject to d20_7:
	0 >= (x20-x7) * (x20-x7) + (y20-y7) * (y20-y7) - 1.0;
subject to d20_8:
	0 >= (x20-x8) * (x20-x8) + (y20-y8) * (y20-y8) - 1.0;
subject to d20_9:
	0 >= (x20-x9) * (x20-x9) + (y20-y9) * (y20-y9) - 1.0;
subject to d20_10:
	0 >= (x20-x10) * (x20-x10) + (y20-y10) * (y20-y10) - 1.0;
subject to d20_11:
	0 >= (x20-x11) * (x20-x11) + (y20-y11) * (y20-y11) - 1.0;
subject to d20_12:
	0 >= (x20-x12) * (x20-x12) + (y20-y12) * (y20-y12) - 1.0;
subject to d20_13:
	0 >= (x20-x13) * (x20-x13) + (y20-y13) * (y20-y13) - 1.0;
subject to d20_14:
	0 >= (x20-x14) * (x20-x14) + (y20-y14) * (y20-y14) - 1.0;
subject to d20_15:
	0 >= (x20-x15) * (x20-x15) + (y20-y15) * (y20-y15) - 1.0;
subject to d20_16:
	0 >= (x20-x16) * (x20-x16) + (y20-y16) * (y20-y16) - 1.0;
subject to d20_17:
	0 >= (x20-x17) * (x20-x17) + (y20-y17) * (y20-y17) - 1.0;
subject to d20_18:
	0 >= (x20-x18) * (x20-x18) + (y20-y18) * (y20-y18) - 1.0;
subject to d20_19:
	0 >= (x20-x19) * (x20-x19) + (y20-y19) * (y20-y19) - 1.0;
subject to d21_0:
	0 >= (x21-x0) * (x21-x0) + (y21-y0) * (y21-y0) - 1.0;
subject to d21_1:
	0 >= (x21-x1) * (x21-x1) + (y21-y1) * (y21-y1) - 1.0;
subject to d21_2:
	0 >= (x21-x2) * (x21-x2) + (y21-y2) * (y21-y2) - 1.0;
subject to d21_3:
	0 >= (x21-x3) * (x21-x3) + (y21-y3) * (y21-y3) - 1.0;
subject to d21_4:
	0 >= (x21-x4) * (x21-x4) + (y21-y4) * (y21-y4) - 1.0;
subject to d21_5:
	0 >= (x21-x5) * (x21-x5) + (y21-y5) * (y21-y5) - 1.0;
subject to d21_6:
	0 >= (x21-x6) * (x21-x6) + (y21-y6) * (y21-y6) - 1.0;
subject to d21_7:
	0 >= (x21-x7) * (x21-x7) + (y21-y7) * (y21-y7) - 1.0;
subject to d21_8:
	0 >= (x21-x8) * (x21-x8) + (y21-y8) * (y21-y8) - 1.0;
subject to d21_9:
	0 >= (x21-x9) * (x21-x9) + (y21-y9) * (y21-y9) - 1.0;
subject to d21_10:
	0 >= (x21-x10) * (x21-x10) + (y21-y10) * (y21-y10) - 1.0;
subject to d21_11:
	0 >= (x21-x11) * (x21-x11) + (y21-y11) * (y21-y11) - 1.0;
subject to d21_12:
	0 >= (x21-x12) * (x21-x12) + (y21-y12) * (y21-y12) - 1.0;
subject to d21_13:
	0 >= (x21-x13) * (x21-x13) + (y21-y13) * (y21-y13) - 1.0;
subject to d21_14:
	0 >= (x21-x14) * (x21-x14) + (y21-y14) * (y21-y14) - 1.0;
subject to d21_15:
	0 >= (x21-x15) * (x21-x15) + (y21-y15) * (y21-y15) - 1.0;
subject to d21_16:
	0 >= (x21-x16) * (x21-x16) + (y21-y16) * (y21-y16) - 1.0;
subject to d21_17:
	0 >= (x21-x17) * (x21-x17) + (y21-y17) * (y21-y17) - 1.0;
subject to d21_18:
	0 >= (x21-x18) * (x21-x18) + (y21-y18) * (y21-y18) - 1.0;
subject to d21_19:
	0 >= (x21-x19) * (x21-x19) + (y21-y19) * (y21-y19) - 1.0;
subject to d21_20:
	0 >= (x21-x20) * (x21-x20) + (y21-y20) * (y21-y20) - 1.0;
subject to d22_0:
	0 >= (x22-x0) * (x22-x0) + (y22-y0) * (y22-y0) - 1.0;
subject to d22_1:
	0 >= (x22-x1) * (x22-x1) + (y22-y1) * (y22-y1) - 1.0;
subject to d22_2:
	0 >= (x22-x2) * (x22-x2) + (y22-y2) * (y22-y2) - 1.0;
subject to d22_3:
	0 >= (x22-x3) * (x22-x3) + (y22-y3) * (y22-y3) - 1.0;
subject to d22_4:
	0 >= (x22-x4) * (x22-x4) + (y22-y4) * (y22-y4) - 1.0;
subject to d22_5:
	0 >= (x22-x5) * (x22-x5) + (y22-y5) * (y22-y5) - 1.0;
subject to d22_6:
	0 >= (x22-x6) * (x22-x6) + (y22-y6) * (y22-y6) - 1.0;
subject to d22_7:
	0 >= (x22-x7) * (x22-x7) + (y22-y7) * (y22-y7) - 1.0;
subject to d22_8:
	0 >= (x22-x8) * (x22-x8) + (y22-y8) * (y22-y8) - 1.0;
subject to d22_9:
	0 >= (x22-x9) * (x22-x9) + (y22-y9) * (y22-y9) - 1.0;
subject to d22_10:
	0 >= (x22-x10) * (x22-x10) + (y22-y10) * (y22-y10) - 1.0;
subject to d22_11:
	0 >= (x22-x11) * (x22-x11) + (y22-y11) * (y22-y11) - 1.0;
subject to d22_12:
	0 >= (x22-x12) * (x22-x12) + (y22-y12) * (y22-y12) - 1.0;
subject to d22_13:
	0 >= (x22-x13) * (x22-x13) + (y22-y13) * (y22-y13) - 1.0;
subject to d22_14:
	0 >= (x22-x14) * (x22-x14) + (y22-y14) * (y22-y14) - 1.0;
subject to d22_15:
	0 >= (x22-x15) * (x22-x15) + (y22-y15) * (y22-y15) - 1.0;
subject to d22_16:
	0 >= (x22-x16) * (x22-x16) + (y22-y16) * (y22-y16) - 1.0;
subject to d22_17:
	0 >= (x22-x17) * (x22-x17) + (y22-y17) * (y22-y17) - 1.0;
subject to d22_18:
	0 >= (x22-x18) * (x22-x18) + (y22-y18) * (y22-y18) - 1.0;
subject to d22_19:
	0 >= (x22-x19) * (x22-x19) + (y22-y19) * (y22-y19) - 1.0;
subject to d22_20:
	0 >= (x22-x20) * (x22-x20) + (y22-y20) * (y22-y20) - 1.0;
subject to d22_21:
	0 >= (x22-x21) * (x22-x21) + (y22-y21) * (y22-y21) - 1.0;
subject to d23_0:
	0 >= (x23-x0) * (x23-x0) + (y23-y0) * (y23-y0) - 1.0;
subject to d23_1:
	0 >= (x23-x1) * (x23-x1) + (y23-y1) * (y23-y1) - 1.0;
subject to d23_2:
	0 >= (x23-x2) * (x23-x2) + (y23-y2) * (y23-y2) - 1.0;
subject to d23_3:
	0 >= (x23-x3) * (x23-x3) + (y23-y3) * (y23-y3) - 1.0;
subject to d23_4:
	0 >= (x23-x4) * (x23-x4) + (y23-y4) * (y23-y4) - 1.0;
subject to d23_5:
	0 >= (x23-x5) * (x23-x5) + (y23-y5) * (y23-y5) - 1.0;
subject to d23_6:
	0 >= (x23-x6) * (x23-x6) + (y23-y6) * (y23-y6) - 1.0;
subject to d23_7:
	0 >= (x23-x7) * (x23-x7) + (y23-y7) * (y23-y7) - 1.0;
subject to d23_8:
	0 >= (x23-x8) * (x23-x8) + (y23-y8) * (y23-y8) - 1.0;
subject to d23_9:
	0 >= (x23-x9) * (x23-x9) + (y23-y9) * (y23-y9) - 1.0;
subject to d23_10:
	0 >= (x23-x10) * (x23-x10) + (y23-y10) * (y23-y10) - 1.0;
subject to d23_11:
	0 >= (x23-x11) * (x23-x11) + (y23-y11) * (y23-y11) - 1.0;
subject to d23_12:
	0 >= (x23-x12) * (x23-x12) + (y23-y12) * (y23-y12) - 1.0;
subject to d23_13:
	0 >= (x23-x13) * (x23-x13) + (y23-y13) * (y23-y13) - 1.0;
subject to d23_14:
	0 >= (x23-x14) * (x23-x14) + (y23-y14) * (y23-y14) - 1.0;
subject to d23_15:
	0 >= (x23-x15) * (x23-x15) + (y23-y15) * (y23-y15) - 1.0;
subject to d23_16:
	0 >= (x23-x16) * (x23-x16) + (y23-y16) * (y23-y16) - 1.0;
subject to d23_17:
	0 >= (x23-x17) * (x23-x17) + (y23-y17) * (y23-y17) - 1.0;
subject to d23_18:
	0 >= (x23-x18) * (x23-x18) + (y23-y18) * (y23-y18) - 1.0;
subject to d23_19:
	0 >= (x23-x19) * (x23-x19) + (y23-y19) * (y23-y19) - 1.0;
subject to d23_20:
	0 >= (x23-x20) * (x23-x20) + (y23-y20) * (y23-y20) - 1.0;
subject to d23_21:
	0 >= (x23-x21) * (x23-x21) + (y23-y21) * (y23-y21) - 1.0;
subject to d23_22:
	0 >= (x23-x22) * (x23-x22) + (y23-y22) * (y23-y22) - 1.0;
subject to d24_0:
	0 >= (x24-x0) * (x24-x0) + (y24-y0) * (y24-y0) - 1.0;
subject to d24_1:
	0 >= (x24-x1) * (x24-x1) + (y24-y1) * (y24-y1) - 1.0;
subject to d24_2:
	0 >= (x24-x2) * (x24-x2) + (y24-y2) * (y24-y2) - 1.0;
subject to d24_3:
	0 >= (x24-x3) * (x24-x3) + (y24-y3) * (y24-y3) - 1.0;
subject to d24_4:
	0 >= (x24-x4) * (x24-x4) + (y24-y4) * (y24-y4) - 1.0;
subject to d24_5:
	0 >= (x24-x5) * (x24-x5) + (y24-y5) * (y24-y5) - 1.0;
subject to d24_6:
	0 >= (x24-x6) * (x24-x6) + (y24-y6) * (y24-y6) - 1.0;
subject to d24_7:
	0 >= (x24-x7) * (x24-x7) + (y24-y7) * (y24-y7) - 1.0;
subject to d24_8:
	0 >= (x24-x8) * (x24-x8) + (y24-y8) * (y24-y8) - 1.0;
subject to d24_9:
	0 >= (x24-x9) * (x24-x9) + (y24-y9) * (y24-y9) - 1.0;
subject to d24_10:
	0 >= (x24-x10) * (x24-x10) + (y24-y10) * (y24-y10) - 1.0;
subject to d24_11:
	0 >= (x24-x11) * (x24-x11) + (y24-y11) * (y24-y11) - 1.0;
subject to d24_12:
	0 >= (x24-x12) * (x24-x12) + (y24-y12) * (y24-y12) - 1.0;
subject to d24_13:
	0 >= (x24-x13) * (x24-x13) + (y24-y13) * (y24-y13) - 1.0;
subject to d24_14:
	0 >= (x24-x14) * (x24-x14) + (y24-y14) * (y24-y14) - 1.0;
subject to d24_15:
	0 >= (x24-x15) * (x24-x15) + (y24-y15) * (y24-y15) - 1.0;
subject to d24_16:
	0 >= (x24-x16) * (x24-x16) + (y24-y16) * (y24-y16) - 1.0;
subject to d24_17:
	0 >= (x24-x17) * (x24-x17) + (y24-y17) * (y24-y17) - 1.0;
subject to d24_18:
	0 >= (x24-x18) * (x24-x18) + (y24-y18) * (y24-y18) - 1.0;
subject to d24_19:
	0 >= (x24-x19) * (x24-x19) + (y24-y19) * (y24-y19) - 1.0;
subject to d24_20:
	0 >= (x24-x20) * (x24-x20) + (y24-y20) * (y24-y20) - 1.0;
subject to d24_21:
	0 >= (x24-x21) * (x24-x21) + (y24-y21) * (y24-y21) - 1.0;
subject to d24_22:
	0 >= (x24-x22) * (x24-x22) + (y24-y22) * (y24-y22) - 1.0;
subject to d24_23:
	0 >= (x24-x23) * (x24-x23) + (y24-y23) * (y24-y23) - 1.0;
subject to d25_0:
	0 >= (x25-x0) * (x25-x0) + (y25-y0) * (y25-y0) - 1.0;
subject to d25_1:
	0 >= (x25-x1) * (x25-x1) + (y25-y1) * (y25-y1) - 1.0;
subject to d25_2:
	0 >= (x25-x2) * (x25-x2) + (y25-y2) * (y25-y2) - 1.0;
subject to d25_3:
	0 >= (x25-x3) * (x25-x3) + (y25-y3) * (y25-y3) - 1.0;
subject to d25_4:
	0 >= (x25-x4) * (x25-x4) + (y25-y4) * (y25-y4) - 1.0;
subject to d25_5:
	0 >= (x25-x5) * (x25-x5) + (y25-y5) * (y25-y5) - 1.0;
subject to d25_6:
	0 >= (x25-x6) * (x25-x6) + (y25-y6) * (y25-y6) - 1.0;
subject to d25_7:
	0 >= (x25-x7) * (x25-x7) + (y25-y7) * (y25-y7) - 1.0;
subject to d25_8:
	0 >= (x25-x8) * (x25-x8) + (y25-y8) * (y25-y8) - 1.0;
subject to d25_9:
	0 >= (x25-x9) * (x25-x9) + (y25-y9) * (y25-y9) - 1.0;
subject to d25_10:
	0 >= (x25-x10) * (x25-x10) + (y25-y10) * (y25-y10) - 1.0;
subject to d25_11:
	0 >= (x25-x11) * (x25-x11) + (y25-y11) * (y25-y11) - 1.0;
subject to d25_12:
	0 >= (x25-x12) * (x25-x12) + (y25-y12) * (y25-y12) - 1.0;
subject to d25_13:
	0 >= (x25-x13) * (x25-x13) + (y25-y13) * (y25-y13) - 1.0;
subject to d25_14:
	0 >= (x25-x14) * (x25-x14) + (y25-y14) * (y25-y14) - 1.0;
subject to d25_15:
	0 >= (x25-x15) * (x25-x15) + (y25-y15) * (y25-y15) - 1.0;
subject to d25_16:
	0 >= (x25-x16) * (x25-x16) + (y25-y16) * (y25-y16) - 1.0;
subject to d25_17:
	0 >= (x25-x17) * (x25-x17) + (y25-y17) * (y25-y17) - 1.0;
subject to d25_18:
	0 >= (x25-x18) * (x25-x18) + (y25-y18) * (y25-y18) - 1.0;
subject to d25_19:
	0 >= (x25-x19) * (x25-x19) + (y25-y19) * (y25-y19) - 1.0;
subject to d25_20:
	0 >= (x25-x20) * (x25-x20) + (y25-y20) * (y25-y20) - 1.0;
subject to d25_21:
	0 >= (x25-x21) * (x25-x21) + (y25-y21) * (y25-y21) - 1.0;
subject to d25_22:
	0 >= (x25-x22) * (x25-x22) + (y25-y22) * (y25-y22) - 1.0;
subject to d25_23:
	0 >= (x25-x23) * (x25-x23) + (y25-y23) * (y25-y23) - 1.0;
subject to d25_24:
	0 >= (x25-x24) * (x25-x24) + (y25-y24) * (y25-y24) - 1.0;
subject to d26_0:
	0 >= (x26-x0) * (x26-x0) + (y26-y0) * (y26-y0) - 1.0;
subject to d26_1:
	0 >= (x26-x1) * (x26-x1) + (y26-y1) * (y26-y1) - 1.0;
subject to d26_2:
	0 >= (x26-x2) * (x26-x2) + (y26-y2) * (y26-y2) - 1.0;
subject to d26_3:
	0 >= (x26-x3) * (x26-x3) + (y26-y3) * (y26-y3) - 1.0;
subject to d26_4:
	0 >= (x26-x4) * (x26-x4) + (y26-y4) * (y26-y4) - 1.0;
subject to d26_5:
	0 >= (x26-x5) * (x26-x5) + (y26-y5) * (y26-y5) - 1.0;
subject to d26_6:
	0 >= (x26-x6) * (x26-x6) + (y26-y6) * (y26-y6) - 1.0;
subject to d26_7:
	0 >= (x26-x7) * (x26-x7) + (y26-y7) * (y26-y7) - 1.0;
subject to d26_8:
	0 >= (x26-x8) * (x26-x8) + (y26-y8) * (y26-y8) - 1.0;
subject to d26_9:
	0 >= (x26-x9) * (x26-x9) + (y26-y9) * (y26-y9) - 1.0;
subject to d26_10:
	0 >= (x26-x10) * (x26-x10) + (y26-y10) * (y26-y10) - 1.0;
subject to d26_11:
	0 >= (x26-x11) * (x26-x11) + (y26-y11) * (y26-y11) - 1.0;
subject to d26_12:
	0 >= (x26-x12) * (x26-x12) + (y26-y12) * (y26-y12) - 1.0;
subject to d26_13:
	0 >= (x26-x13) * (x26-x13) + (y26-y13) * (y26-y13) - 1.0;
subject to d26_14:
	0 >= (x26-x14) * (x26-x14) + (y26-y14) * (y26-y14) - 1.0;
subject to d26_15:
	0 >= (x26-x15) * (x26-x15) + (y26-y15) * (y26-y15) - 1.0;
subject to d26_16:
	0 >= (x26-x16) * (x26-x16) + (y26-y16) * (y26-y16) - 1.0;
subject to d26_17:
	0 >= (x26-x17) * (x26-x17) + (y26-y17) * (y26-y17) - 1.0;
subject to d26_18:
	0 >= (x26-x18) * (x26-x18) + (y26-y18) * (y26-y18) - 1.0;
subject to d26_19:
	0 >= (x26-x19) * (x26-x19) + (y26-y19) * (y26-y19) - 1.0;
subject to d26_20:
	0 >= (x26-x20) * (x26-x20) + (y26-y20) * (y26-y20) - 1.0;
subject to d26_21:
	0 >= (x26-x21) * (x26-x21) + (y26-y21) * (y26-y21) - 1.0;
subject to d26_22:
	0 >= (x26-x22) * (x26-x22) + (y26-y22) * (y26-y22) - 1.0;
subject to d26_23:
	0 >= (x26-x23) * (x26-x23) + (y26-y23) * (y26-y23) - 1.0;
subject to d26_24:
	0 >= (x26-x24) * (x26-x24) + (y26-y24) * (y26-y24) - 1.0;
subject to d26_25:
	0 >= (x26-x25) * (x26-x25) + (y26-y25) * (y26-y25) - 1.0;
subject to d27_0:
	0 >= (x27-x0) * (x27-x0) + (y27-y0) * (y27-y0) - 1.0;
subject to d27_1:
	0 >= (x27-x1) * (x27-x1) + (y27-y1) * (y27-y1) - 1.0;
subject to d27_2:
	0 >= (x27-x2) * (x27-x2) + (y27-y2) * (y27-y2) - 1.0;
subject to d27_3:
	0 >= (x27-x3) * (x27-x3) + (y27-y3) * (y27-y3) - 1.0;
subject to d27_4:
	0 >= (x27-x4) * (x27-x4) + (y27-y4) * (y27-y4) - 1.0;
subject to d27_5:
	0 >= (x27-x5) * (x27-x5) + (y27-y5) * (y27-y5) - 1.0;
subject to d27_6:
	0 >= (x27-x6) * (x27-x6) + (y27-y6) * (y27-y6) - 1.0;
subject to d27_7:
	0 >= (x27-x7) * (x27-x7) + (y27-y7) * (y27-y7) - 1.0;
subject to d27_8:
	0 >= (x27-x8) * (x27-x8) + (y27-y8) * (y27-y8) - 1.0;
subject to d27_9:
	0 >= (x27-x9) * (x27-x9) + (y27-y9) * (y27-y9) - 1.0;
subject to d27_10:
	0 >= (x27-x10) * (x27-x10) + (y27-y10) * (y27-y10) - 1.0;
subject to d27_11:
	0 >= (x27-x11) * (x27-x11) + (y27-y11) * (y27-y11) - 1.0;
subject to d27_12:
	0 >= (x27-x12) * (x27-x12) + (y27-y12) * (y27-y12) - 1.0;
subject to d27_13:
	0 >= (x27-x13) * (x27-x13) + (y27-y13) * (y27-y13) - 1.0;
subject to d27_14:
	0 >= (x27-x14) * (x27-x14) + (y27-y14) * (y27-y14) - 1.0;
subject to d27_15:
	0 >= (x27-x15) * (x27-x15) + (y27-y15) * (y27-y15) - 1.0;
subject to d27_16:
	0 >= (x27-x16) * (x27-x16) + (y27-y16) * (y27-y16) - 1.0;
subject to d27_17:
	0 >= (x27-x17) * (x27-x17) + (y27-y17) * (y27-y17) - 1.0;
subject to d27_18:
	0 >= (x27-x18) * (x27-x18) + (y27-y18) * (y27-y18) - 1.0;
subject to d27_19:
	0 >= (x27-x19) * (x27-x19) + (y27-y19) * (y27-y19) - 1.0;
subject to d27_20:
	0 >= (x27-x20) * (x27-x20) + (y27-y20) * (y27-y20) - 1.0;
subject to d27_21:
	0 >= (x27-x21) * (x27-x21) + (y27-y21) * (y27-y21) - 1.0;
subject to d27_22:
	0 >= (x27-x22) * (x27-x22) + (y27-y22) * (y27-y22) - 1.0;
subject to d27_23:
	0 >= (x27-x23) * (x27-x23) + (y27-y23) * (y27-y23) - 1.0;
subject to d27_24:
	0 >= (x27-x24) * (x27-x24) + (y27-y24) * (y27-y24) - 1.0;
subject to d27_25:
	0 >= (x27-x25) * (x27-x25) + (y27-y25) * (y27-y25) - 1.0;
subject to d27_26:
	0 >= (x27-x26) * (x27-x26) + (y27-y26) * (y27-y26) - 1.0;
subject to d28_0:
	0 >= (x28-x0) * (x28-x0) + (y28-y0) * (y28-y0) - 1.0;
subject to d28_1:
	0 >= (x28-x1) * (x28-x1) + (y28-y1) * (y28-y1) - 1.0;
subject to d28_2:
	0 >= (x28-x2) * (x28-x2) + (y28-y2) * (y28-y2) - 1.0;
subject to d28_3:
	0 >= (x28-x3) * (x28-x3) + (y28-y3) * (y28-y3) - 1.0;
subject to d28_4:
	0 >= (x28-x4) * (x28-x4) + (y28-y4) * (y28-y4) - 1.0;
subject to d28_5:
	0 >= (x28-x5) * (x28-x5) + (y28-y5) * (y28-y5) - 1.0;
subject to d28_6:
	0 >= (x28-x6) * (x28-x6) + (y28-y6) * (y28-y6) - 1.0;
subject to d28_7:
	0 >= (x28-x7) * (x28-x7) + (y28-y7) * (y28-y7) - 1.0;
subject to d28_8:
	0 >= (x28-x8) * (x28-x8) + (y28-y8) * (y28-y8) - 1.0;
subject to d28_9:
	0 >= (x28-x9) * (x28-x9) + (y28-y9) * (y28-y9) - 1.0;
subject to d28_10:
	0 >= (x28-x10) * (x28-x10) + (y28-y10) * (y28-y10) - 1.0;
subject to d28_11:
	0 >= (x28-x11) * (x28-x11) + (y28-y11) * (y28-y11) - 1.0;
subject to d28_12:
	0 >= (x28-x12) * (x28-x12) + (y28-y12) * (y28-y12) - 1.0;
subject to d28_13:
	0 >= (x28-x13) * (x28-x13) + (y28-y13) * (y28-y13) - 1.0;
subject to d28_14:
	0 >= (x28-x14) * (x28-x14) + (y28-y14) * (y28-y14) - 1.0;
subject to d28_15:
	0 >= (x28-x15) * (x28-x15) + (y28-y15) * (y28-y15) - 1.0;
subject to d28_16:
	0 >= (x28-x16) * (x28-x16) + (y28-y16) * (y28-y16) - 1.0;
subject to d28_17:
	0 >= (x28-x17) * (x28-x17) + (y28-y17) * (y28-y17) - 1.0;
subject to d28_18:
	0 >= (x28-x18) * (x28-x18) + (y28-y18) * (y28-y18) - 1.0;
subject to d28_19:
	0 >= (x28-x19) * (x28-x19) + (y28-y19) * (y28-y19) - 1.0;
subject to d28_20:
	0 >= (x28-x20) * (x28-x20) + (y28-y20) * (y28-y20) - 1.0;
subject to d28_21:
	0 >= (x28-x21) * (x28-x21) + (y28-y21) * (y28-y21) - 1.0;
subject to d28_22:
	0 >= (x28-x22) * (x28-x22) + (y28-y22) * (y28-y22) - 1.0;
subject to d28_23:
	0 >= (x28-x23) * (x28-x23) + (y28-y23) * (y28-y23) - 1.0;
subject to d28_24:
	0 >= (x28-x24) * (x28-x24) + (y28-y24) * (y28-y24) - 1.0;
subject to d28_25:
	0 >= (x28-x25) * (x28-x25) + (y28-y25) * (y28-y25) - 1.0;
subject to d28_26:
	0 >= (x28-x26) * (x28-x26) + (y28-y26) * (y28-y26) - 1.0;
subject to d28_27:
	0 >= (x28-x27) * (x28-x27) + (y28-y27) * (y28-y27) - 1.0;
subject to d29_0:
	0 >= (x29-x0) * (x29-x0) + (y29-y0) * (y29-y0) - 1.0;
subject to d29_1:
	0 >= (x29-x1) * (x29-x1) + (y29-y1) * (y29-y1) - 1.0;
subject to d29_2:
	0 >= (x29-x2) * (x29-x2) + (y29-y2) * (y29-y2) - 1.0;
subject to d29_3:
	0 >= (x29-x3) * (x29-x3) + (y29-y3) * (y29-y3) - 1.0;
subject to d29_4:
	0 >= (x29-x4) * (x29-x4) + (y29-y4) * (y29-y4) - 1.0;
subject to d29_5:
	0 >= (x29-x5) * (x29-x5) + (y29-y5) * (y29-y5) - 1.0;
subject to d29_6:
	0 >= (x29-x6) * (x29-x6) + (y29-y6) * (y29-y6) - 1.0;
subject to d29_7:
	0 >= (x29-x7) * (x29-x7) + (y29-y7) * (y29-y7) - 1.0;
subject to d29_8:
	0 >= (x29-x8) * (x29-x8) + (y29-y8) * (y29-y8) - 1.0;
subject to d29_9:
	0 >= (x29-x9) * (x29-x9) + (y29-y9) * (y29-y9) - 1.0;
subject to d29_10:
	0 >= (x29-x10) * (x29-x10) + (y29-y10) * (y29-y10) - 1.0;
subject to d29_11:
	0 >= (x29-x11) * (x29-x11) + (y29-y11) * (y29-y11) - 1.0;
subject to d29_12:
	0 >= (x29-x12) * (x29-x12) + (y29-y12) * (y29-y12) - 1.0;
subject to d29_13:
	0 >= (x29-x13) * (x29-x13) + (y29-y13) * (y29-y13) - 1.0;
subject to d29_14:
	0 >= (x29-x14) * (x29-x14) + (y29-y14) * (y29-y14) - 1.0;
subject to d29_15:
	0 >= (x29-x15) * (x29-x15) + (y29-y15) * (y29-y15) - 1.0;
subject to d29_16:
	0 >= (x29-x16) * (x29-x16) + (y29-y16) * (y29-y16) - 1.0;
subject to d29_17:
	0 >= (x29-x17) * (x29-x17) + (y29-y17) * (y29-y17) - 1.0;
subject to d29_18:
	0 >= (x29-x18) * (x29-x18) + (y29-y18) * (y29-y18) - 1.0;
subject to d29_19:
	0 >= (x29-x19) * (x29-x19) + (y29-y19) * (y29-y19) - 1.0;
subject to d29_20:
	0 >= (x29-x20) * (x29-x20) + (y29-y20) * (y29-y20) - 1.0;
subject to d29_21:
	0 >= (x29-x21) * (x29-x21) + (y29-y21) * (y29-y21) - 1.0;
subject to d29_22:
	0 >= (x29-x22) * (x29-x22) + (y29-y22) * (y29-y22) - 1.0;
subject to d29_23:
	0 >= (x29-x23) * (x29-x23) + (y29-y23) * (y29-y23) - 1.0;
subject to d29_24:
	0 >= (x29-x24) * (x29-x24) + (y29-y24) * (y29-y24) - 1.0;
subject to d29_25:
	0 >= (x29-x25) * (x29-x25) + (y29-y25) * (y29-y25) - 1.0;
subject to d29_26:
	0 >= (x29-x26) * (x29-x26) + (y29-y26) * (y29-y26) - 1.0;
subject to d29_27:
	0 >= (x29-x27) * (x29-x27) + (y29-y27) * (y29-y27) - 1.0;
subject to d29_28:
	0 >= (x29-x28) * (x29-x28) + (y29-y28) * (y29-y28) - 1.0;
subject to d30_0:
	0 >= (x30-x0) * (x30-x0) + (y30-y0) * (y30-y0) - 1.0;
subject to d30_1:
	0 >= (x30-x1) * (x30-x1) + (y30-y1) * (y30-y1) - 1.0;
subject to d30_2:
	0 >= (x30-x2) * (x30-x2) + (y30-y2) * (y30-y2) - 1.0;
subject to d30_3:
	0 >= (x30-x3) * (x30-x3) + (y30-y3) * (y30-y3) - 1.0;
subject to d30_4:
	0 >= (x30-x4) * (x30-x4) + (y30-y4) * (y30-y4) - 1.0;
subject to d30_5:
	0 >= (x30-x5) * (x30-x5) + (y30-y5) * (y30-y5) - 1.0;
subject to d30_6:
	0 >= (x30-x6) * (x30-x6) + (y30-y6) * (y30-y6) - 1.0;
subject to d30_7:
	0 >= (x30-x7) * (x30-x7) + (y30-y7) * (y30-y7) - 1.0;
subject to d30_8:
	0 >= (x30-x8) * (x30-x8) + (y30-y8) * (y30-y8) - 1.0;
subject to d30_9:
	0 >= (x30-x9) * (x30-x9) + (y30-y9) * (y30-y9) - 1.0;
subject to d30_10:
	0 >= (x30-x10) * (x30-x10) + (y30-y10) * (y30-y10) - 1.0;
subject to d30_11:
	0 >= (x30-x11) * (x30-x11) + (y30-y11) * (y30-y11) - 1.0;
subject to d30_12:
	0 >= (x30-x12) * (x30-x12) + (y30-y12) * (y30-y12) - 1.0;
subject to d30_13:
	0 >= (x30-x13) * (x30-x13) + (y30-y13) * (y30-y13) - 1.0;
subject to d30_14:
	0 >= (x30-x14) * (x30-x14) + (y30-y14) * (y30-y14) - 1.0;
subject to d30_15:
	0 >= (x30-x15) * (x30-x15) + (y30-y15) * (y30-y15) - 1.0;
subject to d30_16:
	0 >= (x30-x16) * (x30-x16) + (y30-y16) * (y30-y16) - 1.0;
subject to d30_17:
	0 >= (x30-x17) * (x30-x17) + (y30-y17) * (y30-y17) - 1.0;
subject to d30_18:
	0 >= (x30-x18) * (x30-x18) + (y30-y18) * (y30-y18) - 1.0;
subject to d30_19:
	0 >= (x30-x19) * (x30-x19) + (y30-y19) * (y30-y19) - 1.0;
subject to d30_20:
	0 >= (x30-x20) * (x30-x20) + (y30-y20) * (y30-y20) - 1.0;
subject to d30_21:
	0 >= (x30-x21) * (x30-x21) + (y30-y21) * (y30-y21) - 1.0;
subject to d30_22:
	0 >= (x30-x22) * (x30-x22) + (y30-y22) * (y30-y22) - 1.0;
subject to d30_23:
	0 >= (x30-x23) * (x30-x23) + (y30-y23) * (y30-y23) - 1.0;
subject to d30_24:
	0 >= (x30-x24) * (x30-x24) + (y30-y24) * (y30-y24) - 1.0;
subject to d30_25:
	0 >= (x30-x25) * (x30-x25) + (y30-y25) * (y30-y25) - 1.0;
subject to d30_26:
	0 >= (x30-x26) * (x30-x26) + (y30-y26) * (y30-y26) - 1.0;
subject to d30_27:
	0 >= (x30-x27) * (x30-x27) + (y30-y27) * (y30-y27) - 1.0;
subject to d30_28:
	0 >= (x30-x28) * (x30-x28) + (y30-y28) * (y30-y28) - 1.0;
subject to d30_29:
	0 >= (x30-x29) * (x30-x29) + (y30-y29) * (y30-y29) - 1.0;
subject to d31_0:
	0 >= (x31-x0) * (x31-x0) + (y31-y0) * (y31-y0) - 1.0;
subject to d31_1:
	0 >= (x31-x1) * (x31-x1) + (y31-y1) * (y31-y1) - 1.0;
subject to d31_2:
	0 >= (x31-x2) * (x31-x2) + (y31-y2) * (y31-y2) - 1.0;
subject to d31_3:
	0 >= (x31-x3) * (x31-x3) + (y31-y3) * (y31-y3) - 1.0;
subject to d31_4:
	0 >= (x31-x4) * (x31-x4) + (y31-y4) * (y31-y4) - 1.0;
subject to d31_5:
	0 >= (x31-x5) * (x31-x5) + (y31-y5) * (y31-y5) - 1.0;
subject to d31_6:
	0 >= (x31-x6) * (x31-x6) + (y31-y6) * (y31-y6) - 1.0;
subject to d31_7:
	0 >= (x31-x7) * (x31-x7) + (y31-y7) * (y31-y7) - 1.0;
subject to d31_8:
	0 >= (x31-x8) * (x31-x8) + (y31-y8) * (y31-y8) - 1.0;
subject to d31_9:
	0 >= (x31-x9) * (x31-x9) + (y31-y9) * (y31-y9) - 1.0;
subject to d31_10:
	0 >= (x31-x10) * (x31-x10) + (y31-y10) * (y31-y10) - 1.0;
subject to d31_11:
	0 >= (x31-x11) * (x31-x11) + (y31-y11) * (y31-y11) - 1.0;
subject to d31_12:
	0 >= (x31-x12) * (x31-x12) + (y31-y12) * (y31-y12) - 1.0;
subject to d31_13:
	0 >= (x31-x13) * (x31-x13) + (y31-y13) * (y31-y13) - 1.0;
subject to d31_14:
	0 >= (x31-x14) * (x31-x14) + (y31-y14) * (y31-y14) - 1.0;
subject to d31_15:
	0 >= (x31-x15) * (x31-x15) + (y31-y15) * (y31-y15) - 1.0;
subject to d31_16:
	0 >= (x31-x16) * (x31-x16) + (y31-y16) * (y31-y16) - 1.0;
subject to d31_17:
	0 >= (x31-x17) * (x31-x17) + (y31-y17) * (y31-y17) - 1.0;
subject to d31_18:
	0 >= (x31-x18) * (x31-x18) + (y31-y18) * (y31-y18) - 1.0;
subject to d31_19:
	0 >= (x31-x19) * (x31-x19) + (y31-y19) * (y31-y19) - 1.0;
subject to d31_20:
	0 >= (x31-x20) * (x31-x20) + (y31-y20) * (y31-y20) - 1.0;
subject to d31_21:
	0 >= (x31-x21) * (x31-x21) + (y31-y21) * (y31-y21) - 1.0;
subject to d31_22:
	0 >= (x31-x22) * (x31-x22) + (y31-y22) * (y31-y22) - 1.0;
subject to d31_23:
	0 >= (x31-x23) * (x31-x23) + (y31-y23) * (y31-y23) - 1.0;
subject to d31_24:
	0 >= (x31-x24) * (x31-x24) + (y31-y24) * (y31-y24) - 1.0;
subject to d31_25:
	0 >= (x31-x25) * (x31-x25) + (y31-y25) * (y31-y25) - 1.0;
subject to d31_26:
	0 >= (x31-x26) * (x31-x26) + (y31-y26) * (y31-y26) - 1.0;
subject to d31_27:
	0 >= (x31-x27) * (x31-x27) + (y31-y27) * (y31-y27) - 1.0;
subject to d31_28:
	0 >= (x31-x28) * (x31-x28) + (y31-y28) * (y31-y28) - 1.0;
subject to d31_29:
	0 >= (x31-x29) * (x31-x29) + (y31-y29) * (y31-y29) - 1.0;
subject to d31_30:
	0 >= (x31-x30) * (x31-x30) + (y31-y30) * (y31-y30) - 1.0;
subject to d32_0:
	0 >= (x32-x0) * (x32-x0) + (y32-y0) * (y32-y0) - 1.0;
subject to d32_1:
	0 >= (x32-x1) * (x32-x1) + (y32-y1) * (y32-y1) - 1.0;
subject to d32_2:
	0 >= (x32-x2) * (x32-x2) + (y32-y2) * (y32-y2) - 1.0;
subject to d32_3:
	0 >= (x32-x3) * (x32-x3) + (y32-y3) * (y32-y3) - 1.0;
subject to d32_4:
	0 >= (x32-x4) * (x32-x4) + (y32-y4) * (y32-y4) - 1.0;
subject to d32_5:
	0 >= (x32-x5) * (x32-x5) + (y32-y5) * (y32-y5) - 1.0;
subject to d32_6:
	0 >= (x32-x6) * (x32-x6) + (y32-y6) * (y32-y6) - 1.0;
subject to d32_7:
	0 >= (x32-x7) * (x32-x7) + (y32-y7) * (y32-y7) - 1.0;
subject to d32_8:
	0 >= (x32-x8) * (x32-x8) + (y32-y8) * (y32-y8) - 1.0;
subject to d32_9:
	0 >= (x32-x9) * (x32-x9) + (y32-y9) * (y32-y9) - 1.0;
subject to d32_10:
	0 >= (x32-x10) * (x32-x10) + (y32-y10) * (y32-y10) - 1.0;
subject to d32_11:
	0 >= (x32-x11) * (x32-x11) + (y32-y11) * (y32-y11) - 1.0;
subject to d32_12:
	0 >= (x32-x12) * (x32-x12) + (y32-y12) * (y32-y12) - 1.0;
subject to d32_13:
	0 >= (x32-x13) * (x32-x13) + (y32-y13) * (y32-y13) - 1.0;
subject to d32_14:
	0 >= (x32-x14) * (x32-x14) + (y32-y14) * (y32-y14) - 1.0;
subject to d32_15:
	0 >= (x32-x15) * (x32-x15) + (y32-y15) * (y32-y15) - 1.0;
subject to d32_16:
	0 >= (x32-x16) * (x32-x16) + (y32-y16) * (y32-y16) - 1.0;
subject to d32_17:
	0 >= (x32-x17) * (x32-x17) + (y32-y17) * (y32-y17) - 1.0;
subject to d32_18:
	0 >= (x32-x18) * (x32-x18) + (y32-y18) * (y32-y18) - 1.0;
subject to d32_19:
	0 >= (x32-x19) * (x32-x19) + (y32-y19) * (y32-y19) - 1.0;
subject to d32_20:
	0 >= (x32-x20) * (x32-x20) + (y32-y20) * (y32-y20) - 1.0;
subject to d32_21:
	0 >= (x32-x21) * (x32-x21) + (y32-y21) * (y32-y21) - 1.0;
subject to d32_22:
	0 >= (x32-x22) * (x32-x22) + (y32-y22) * (y32-y22) - 1.0;
subject to d32_23:
	0 >= (x32-x23) * (x32-x23) + (y32-y23) * (y32-y23) - 1.0;
subject to d32_24:
	0 >= (x32-x24) * (x32-x24) + (y32-y24) * (y32-y24) - 1.0;
subject to d32_25:
	0 >= (x32-x25) * (x32-x25) + (y32-y25) * (y32-y25) - 1.0;
subject to d32_26:
	0 >= (x32-x26) * (x32-x26) + (y32-y26) * (y32-y26) - 1.0;
subject to d32_27:
	0 >= (x32-x27) * (x32-x27) + (y32-y27) * (y32-y27) - 1.0;
subject to d32_28:
	0 >= (x32-x28) * (x32-x28) + (y32-y28) * (y32-y28) - 1.0;
subject to d32_29:
	0 >= (x32-x29) * (x32-x29) + (y32-y29) * (y32-y29) - 1.0;
subject to d32_30:
	0 >= (x32-x30) * (x32-x30) + (y32-y30) * (y32-y30) - 1.0;
subject to d32_31:
	0 >= (x32-x31) * (x32-x31) + (y32-y31) * (y32-y31) - 1.0;
subject to d33_0:
	0 >= (x33-x0) * (x33-x0) + (y33-y0) * (y33-y0) - 1.0;
subject to d33_1:
	0 >= (x33-x1) * (x33-x1) + (y33-y1) * (y33-y1) - 1.0;
subject to d33_2:
	0 >= (x33-x2) * (x33-x2) + (y33-y2) * (y33-y2) - 1.0;
subject to d33_3:
	0 >= (x33-x3) * (x33-x3) + (y33-y3) * (y33-y3) - 1.0;
subject to d33_4:
	0 >= (x33-x4) * (x33-x4) + (y33-y4) * (y33-y4) - 1.0;
subject to d33_5:
	0 >= (x33-x5) * (x33-x5) + (y33-y5) * (y33-y5) - 1.0;
subject to d33_6:
	0 >= (x33-x6) * (x33-x6) + (y33-y6) * (y33-y6) - 1.0;
subject to d33_7:
	0 >= (x33-x7) * (x33-x7) + (y33-y7) * (y33-y7) - 1.0;
subject to d33_8:
	0 >= (x33-x8) * (x33-x8) + (y33-y8) * (y33-y8) - 1.0;
subject to d33_9:
	0 >= (x33-x9) * (x33-x9) + (y33-y9) * (y33-y9) - 1.0;
subject to d33_10:
	0 >= (x33-x10) * (x33-x10) + (y33-y10) * (y33-y10) - 1.0;
subject to d33_11:
	0 >= (x33-x11) * (x33-x11) + (y33-y11) * (y33-y11) - 1.0;
subject to d33_12:
	0 >= (x33-x12) * (x33-x12) + (y33-y12) * (y33-y12) - 1.0;
subject to d33_13:
	0 >= (x33-x13) * (x33-x13) + (y33-y13) * (y33-y13) - 1.0;
subject to d33_14:
	0 >= (x33-x14) * (x33-x14) + (y33-y14) * (y33-y14) - 1.0;
subject to d33_15:
	0 >= (x33-x15) * (x33-x15) + (y33-y15) * (y33-y15) - 1.0;
subject to d33_16:
	0 >= (x33-x16) * (x33-x16) + (y33-y16) * (y33-y16) - 1.0;
subject to d33_17:
	0 >= (x33-x17) * (x33-x17) + (y33-y17) * (y33-y17) - 1.0;
subject to d33_18:
	0 >= (x33-x18) * (x33-x18) + (y33-y18) * (y33-y18) - 1.0;
subject to d33_19:
	0 >= (x33-x19) * (x33-x19) + (y33-y19) * (y33-y19) - 1.0;
subject to d33_20:
	0 >= (x33-x20) * (x33-x20) + (y33-y20) * (y33-y20) - 1.0;
subject to d33_21:
	0 >= (x33-x21) * (x33-x21) + (y33-y21) * (y33-y21) - 1.0;
subject to d33_22:
	0 >= (x33-x22) * (x33-x22) + (y33-y22) * (y33-y22) - 1.0;
subject to d33_23:
	0 >= (x33-x23) * (x33-x23) + (y33-y23) * (y33-y23) - 1.0;
subject to d33_24:
	0 >= (x33-x24) * (x33-x24) + (y33-y24) * (y33-y24) - 1.0;
subject to d33_25:
	0 >= (x33-x25) * (x33-x25) + (y33-y25) * (y33-y25) - 1.0;
subject to d33_26:
	0 >= (x33-x26) * (x33-x26) + (y33-y26) * (y33-y26) - 1.0;
subject to d33_27:
	0 >= (x33-x27) * (x33-x27) + (y33-y27) * (y33-y27) - 1.0;
subject to d33_28:
	0 >= (x33-x28) * (x33-x28) + (y33-y28) * (y33-y28) - 1.0;
subject to d33_29:
	0 >= (x33-x29) * (x33-x29) + (y33-y29) * (y33-y29) - 1.0;
subject to d33_30:
	0 >= (x33-x30) * (x33-x30) + (y33-y30) * (y33-y30) - 1.0;
subject to d33_31:
	0 >= (x33-x31) * (x33-x31) + (y33-y31) * (y33-y31) - 1.0;
subject to d33_32:
	0 >= (x33-x32) * (x33-x32) + (y33-y32) * (y33-y32) - 1.0;
subject to d34_0:
	0 >= (x34-x0) * (x34-x0) + (y34-y0) * (y34-y0) - 1.0;
subject to d34_1:
	0 >= (x34-x1) * (x34-x1) + (y34-y1) * (y34-y1) - 1.0;
subject to d34_2:
	0 >= (x34-x2) * (x34-x2) + (y34-y2) * (y34-y2) - 1.0;
subject to d34_3:
	0 >= (x34-x3) * (x34-x3) + (y34-y3) * (y34-y3) - 1.0;
subject to d34_4:
	0 >= (x34-x4) * (x34-x4) + (y34-y4) * (y34-y4) - 1.0;
subject to d34_5:
	0 >= (x34-x5) * (x34-x5) + (y34-y5) * (y34-y5) - 1.0;
subject to d34_6:
	0 >= (x34-x6) * (x34-x6) + (y34-y6) * (y34-y6) - 1.0;
subject to d34_7:
	0 >= (x34-x7) * (x34-x7) + (y34-y7) * (y34-y7) - 1.0;
subject to d34_8:
	0 >= (x34-x8) * (x34-x8) + (y34-y8) * (y34-y8) - 1.0;
subject to d34_9:
	0 >= (x34-x9) * (x34-x9) + (y34-y9) * (y34-y9) - 1.0;
subject to d34_10:
	0 >= (x34-x10) * (x34-x10) + (y34-y10) * (y34-y10) - 1.0;
subject to d34_11:
	0 >= (x34-x11) * (x34-x11) + (y34-y11) * (y34-y11) - 1.0;
subject to d34_12:
	0 >= (x34-x12) * (x34-x12) + (y34-y12) * (y34-y12) - 1.0;
subject to d34_13:
	0 >= (x34-x13) * (x34-x13) + (y34-y13) * (y34-y13) - 1.0;
subject to d34_14:
	0 >= (x34-x14) * (x34-x14) + (y34-y14) * (y34-y14) - 1.0;
subject to d34_15:
	0 >= (x34-x15) * (x34-x15) + (y34-y15) * (y34-y15) - 1.0;
subject to d34_16:
	0 >= (x34-x16) * (x34-x16) + (y34-y16) * (y34-y16) - 1.0;
subject to d34_17:
	0 >= (x34-x17) * (x34-x17) + (y34-y17) * (y34-y17) - 1.0;
subject to d34_18:
	0 >= (x34-x18) * (x34-x18) + (y34-y18) * (y34-y18) - 1.0;
subject to d34_19:
	0 >= (x34-x19) * (x34-x19) + (y34-y19) * (y34-y19) - 1.0;
subject to d34_20:
	0 >= (x34-x20) * (x34-x20) + (y34-y20) * (y34-y20) - 1.0;
subject to d34_21:
	0 >= (x34-x21) * (x34-x21) + (y34-y21) * (y34-y21) - 1.0;
subject to d34_22:
	0 >= (x34-x22) * (x34-x22) + (y34-y22) * (y34-y22) - 1.0;
subject to d34_23:
	0 >= (x34-x23) * (x34-x23) + (y34-y23) * (y34-y23) - 1.0;
subject to d34_24:
	0 >= (x34-x24) * (x34-x24) + (y34-y24) * (y34-y24) - 1.0;
subject to d34_25:
	0 >= (x34-x25) * (x34-x25) + (y34-y25) * (y34-y25) - 1.0;
subject to d34_26:
	0 >= (x34-x26) * (x34-x26) + (y34-y26) * (y34-y26) - 1.0;
subject to d34_27:
	0 >= (x34-x27) * (x34-x27) + (y34-y27) * (y34-y27) - 1.0;
subject to d34_28:
	0 >= (x34-x28) * (x34-x28) + (y34-y28) * (y34-y28) - 1.0;
subject to d34_29:
	0 >= (x34-x29) * (x34-x29) + (y34-y29) * (y34-y29) - 1.0;
subject to d34_30:
	0 >= (x34-x30) * (x34-x30) + (y34-y30) * (y34-y30) - 1.0;
subject to d34_31:
	0 >= (x34-x31) * (x34-x31) + (y34-y31) * (y34-y31) - 1.0;
subject to d34_32:
	0 >= (x34-x32) * (x34-x32) + (y34-y32) * (y34-y32) - 1.0;
subject to d34_33:
	0 >= (x34-x33) * (x34-x33) + (y34-y33) * (y34-y33) - 1.0;
subject to d35_0:
	0 >= (x35-x0) * (x35-x0) + (y35-y0) * (y35-y0) - 1.0;
subject to d35_1:
	0 >= (x35-x1) * (x35-x1) + (y35-y1) * (y35-y1) - 1.0;
subject to d35_2:
	0 >= (x35-x2) * (x35-x2) + (y35-y2) * (y35-y2) - 1.0;
subject to d35_3:
	0 >= (x35-x3) * (x35-x3) + (y35-y3) * (y35-y3) - 1.0;
subject to d35_4:
	0 >= (x35-x4) * (x35-x4) + (y35-y4) * (y35-y4) - 1.0;
subject to d35_5:
	0 >= (x35-x5) * (x35-x5) + (y35-y5) * (y35-y5) - 1.0;
subject to d35_6:
	0 >= (x35-x6) * (x35-x6) + (y35-y6) * (y35-y6) - 1.0;
subject to d35_7:
	0 >= (x35-x7) * (x35-x7) + (y35-y7) * (y35-y7) - 1.0;
subject to d35_8:
	0 >= (x35-x8) * (x35-x8) + (y35-y8) * (y35-y8) - 1.0;
subject to d35_9:
	0 >= (x35-x9) * (x35-x9) + (y35-y9) * (y35-y9) - 1.0;
subject to d35_10:
	0 >= (x35-x10) * (x35-x10) + (y35-y10) * (y35-y10) - 1.0;
subject to d35_11:
	0 >= (x35-x11) * (x35-x11) + (y35-y11) * (y35-y11) - 1.0;
subject to d35_12:
	0 >= (x35-x12) * (x35-x12) + (y35-y12) * (y35-y12) - 1.0;
subject to d35_13:
	0 >= (x35-x13) * (x35-x13) + (y35-y13) * (y35-y13) - 1.0;
subject to d35_14:
	0 >= (x35-x14) * (x35-x14) + (y35-y14) * (y35-y14) - 1.0;
subject to d35_15:
	0 >= (x35-x15) * (x35-x15) + (y35-y15) * (y35-y15) - 1.0;
subject to d35_16:
	0 >= (x35-x16) * (x35-x16) + (y35-y16) * (y35-y16) - 1.0;
subject to d35_17:
	0 >= (x35-x17) * (x35-x17) + (y35-y17) * (y35-y17) - 1.0;
subject to d35_18:
	0 >= (x35-x18) * (x35-x18) + (y35-y18) * (y35-y18) - 1.0;
subject to d35_19:
	0 >= (x35-x19) * (x35-x19) + (y35-y19) * (y35-y19) - 1.0;
subject to d35_20:
	0 >= (x35-x20) * (x35-x20) + (y35-y20) * (y35-y20) - 1.0;
subject to d35_21:
	0 >= (x35-x21) * (x35-x21) + (y35-y21) * (y35-y21) - 1.0;
subject to d35_22:
	0 >= (x35-x22) * (x35-x22) + (y35-y22) * (y35-y22) - 1.0;
subject to d35_23:
	0 >= (x35-x23) * (x35-x23) + (y35-y23) * (y35-y23) - 1.0;
subject to d35_24:
	0 >= (x35-x24) * (x35-x24) + (y35-y24) * (y35-y24) - 1.0;
subject to d35_25:
	0 >= (x35-x25) * (x35-x25) + (y35-y25) * (y35-y25) - 1.0;
subject to d35_26:
	0 >= (x35-x26) * (x35-x26) + (y35-y26) * (y35-y26) - 1.0;
subject to d35_27:
	0 >= (x35-x27) * (x35-x27) + (y35-y27) * (y35-y27) - 1.0;
subject to d35_28:
	0 >= (x35-x28) * (x35-x28) + (y35-y28) * (y35-y28) - 1.0;
subject to d35_29:
	0 >= (x35-x29) * (x35-x29) + (y35-y29) * (y35-y29) - 1.0;
subject to d35_30:
	0 >= (x35-x30) * (x35-x30) + (y35-y30) * (y35-y30) - 1.0;
subject to d35_31:
	0 >= (x35-x31) * (x35-x31) + (y35-y31) * (y35-y31) - 1.0;
subject to d35_32:
	0 >= (x35-x32) * (x35-x32) + (y35-y32) * (y35-y32) - 1.0;
subject to d35_33:
	0 >= (x35-x33) * (x35-x33) + (y35-y33) * (y35-y33) - 1.0;
subject to d35_34:
	0 >= (x35-x34) * (x35-x34) + (y35-y34) * (y35-y34) - 1.0;
subject to d36_0:
	0 >= (x36-x0) * (x36-x0) + (y36-y0) * (y36-y0) - 1.0;
subject to d36_1:
	0 >= (x36-x1) * (x36-x1) + (y36-y1) * (y36-y1) - 1.0;
subject to d36_2:
	0 >= (x36-x2) * (x36-x2) + (y36-y2) * (y36-y2) - 1.0;
subject to d36_3:
	0 >= (x36-x3) * (x36-x3) + (y36-y3) * (y36-y3) - 1.0;
subject to d36_4:
	0 >= (x36-x4) * (x36-x4) + (y36-y4) * (y36-y4) - 1.0;
subject to d36_5:
	0 >= (x36-x5) * (x36-x5) + (y36-y5) * (y36-y5) - 1.0;
subject to d36_6:
	0 >= (x36-x6) * (x36-x6) + (y36-y6) * (y36-y6) - 1.0;
subject to d36_7:
	0 >= (x36-x7) * (x36-x7) + (y36-y7) * (y36-y7) - 1.0;
subject to d36_8:
	0 >= (x36-x8) * (x36-x8) + (y36-y8) * (y36-y8) - 1.0;
subject to d36_9:
	0 >= (x36-x9) * (x36-x9) + (y36-y9) * (y36-y9) - 1.0;
subject to d36_10:
	0 >= (x36-x10) * (x36-x10) + (y36-y10) * (y36-y10) - 1.0;
subject to d36_11:
	0 >= (x36-x11) * (x36-x11) + (y36-y11) * (y36-y11) - 1.0;
subject to d36_12:
	0 >= (x36-x12) * (x36-x12) + (y36-y12) * (y36-y12) - 1.0;
subject to d36_13:
	0 >= (x36-x13) * (x36-x13) + (y36-y13) * (y36-y13) - 1.0;
subject to d36_14:
	0 >= (x36-x14) * (x36-x14) + (y36-y14) * (y36-y14) - 1.0;
subject to d36_15:
	0 >= (x36-x15) * (x36-x15) + (y36-y15) * (y36-y15) - 1.0;
subject to d36_16:
	0 >= (x36-x16) * (x36-x16) + (y36-y16) * (y36-y16) - 1.0;
subject to d36_17:
	0 >= (x36-x17) * (x36-x17) + (y36-y17) * (y36-y17) - 1.0;
subject to d36_18:
	0 >= (x36-x18) * (x36-x18) + (y36-y18) * (y36-y18) - 1.0;
subject to d36_19:
	0 >= (x36-x19) * (x36-x19) + (y36-y19) * (y36-y19) - 1.0;
subject to d36_20:
	0 >= (x36-x20) * (x36-x20) + (y36-y20) * (y36-y20) - 1.0;
subject to d36_21:
	0 >= (x36-x21) * (x36-x21) + (y36-y21) * (y36-y21) - 1.0;
subject to d36_22:
	0 >= (x36-x22) * (x36-x22) + (y36-y22) * (y36-y22) - 1.0;
subject to d36_23:
	0 >= (x36-x23) * (x36-x23) + (y36-y23) * (y36-y23) - 1.0;
subject to d36_24:
	0 >= (x36-x24) * (x36-x24) + (y36-y24) * (y36-y24) - 1.0;
subject to d36_25:
	0 >= (x36-x25) * (x36-x25) + (y36-y25) * (y36-y25) - 1.0;
subject to d36_26:
	0 >= (x36-x26) * (x36-x26) + (y36-y26) * (y36-y26) - 1.0;
subject to d36_27:
	0 >= (x36-x27) * (x36-x27) + (y36-y27) * (y36-y27) - 1.0;
subject to d36_28:
	0 >= (x36-x28) * (x36-x28) + (y36-y28) * (y36-y28) - 1.0;
subject to d36_29:
	0 >= (x36-x29) * (x36-x29) + (y36-y29) * (y36-y29) - 1.0;
subject to d36_30:
	0 >= (x36-x30) * (x36-x30) + (y36-y30) * (y36-y30) - 1.0;
subject to d36_31:
	0 >= (x36-x31) * (x36-x31) + (y36-y31) * (y36-y31) - 1.0;
subject to d36_32:
	0 >= (x36-x32) * (x36-x32) + (y36-y32) * (y36-y32) - 1.0;
subject to d36_33:
	0 >= (x36-x33) * (x36-x33) + (y36-y33) * (y36-y33) - 1.0;
subject to d36_34:
	0 >= (x36-x34) * (x36-x34) + (y36-y34) * (y36-y34) - 1.0;
subject to d36_35:
	0 >= (x36-x35) * (x36-x35) + (y36-y35) * (y36-y35) - 1.0;
subject to d37_0:
	0 >= (x37-x0) * (x37-x0) + (y37-y0) * (y37-y0) - 1.0;
subject to d37_1:
	0 >= (x37-x1) * (x37-x1) + (y37-y1) * (y37-y1) - 1.0;
subject to d37_2:
	0 >= (x37-x2) * (x37-x2) + (y37-y2) * (y37-y2) - 1.0;
subject to d37_3:
	0 >= (x37-x3) * (x37-x3) + (y37-y3) * (y37-y3) - 1.0;
subject to d37_4:
	0 >= (x37-x4) * (x37-x4) + (y37-y4) * (y37-y4) - 1.0;
subject to d37_5:
	0 >= (x37-x5) * (x37-x5) + (y37-y5) * (y37-y5) - 1.0;
subject to d37_6:
	0 >= (x37-x6) * (x37-x6) + (y37-y6) * (y37-y6) - 1.0;
subject to d37_7:
	0 >= (x37-x7) * (x37-x7) + (y37-y7) * (y37-y7) - 1.0;
subject to d37_8:
	0 >= (x37-x8) * (x37-x8) + (y37-y8) * (y37-y8) - 1.0;
subject to d37_9:
	0 >= (x37-x9) * (x37-x9) + (y37-y9) * (y37-y9) - 1.0;
subject to d37_10:
	0 >= (x37-x10) * (x37-x10) + (y37-y10) * (y37-y10) - 1.0;
subject to d37_11:
	0 >= (x37-x11) * (x37-x11) + (y37-y11) * (y37-y11) - 1.0;
subject to d37_12:
	0 >= (x37-x12) * (x37-x12) + (y37-y12) * (y37-y12) - 1.0;
subject to d37_13:
	0 >= (x37-x13) * (x37-x13) + (y37-y13) * (y37-y13) - 1.0;
subject to d37_14:
	0 >= (x37-x14) * (x37-x14) + (y37-y14) * (y37-y14) - 1.0;
subject to d37_15:
	0 >= (x37-x15) * (x37-x15) + (y37-y15) * (y37-y15) - 1.0;
subject to d37_16:
	0 >= (x37-x16) * (x37-x16) + (y37-y16) * (y37-y16) - 1.0;
subject to d37_17:
	0 >= (x37-x17) * (x37-x17) + (y37-y17) * (y37-y17) - 1.0;
subject to d37_18:
	0 >= (x37-x18) * (x37-x18) + (y37-y18) * (y37-y18) - 1.0;
subject to d37_19:
	0 >= (x37-x19) * (x37-x19) + (y37-y19) * (y37-y19) - 1.0;
subject to d37_20:
	0 >= (x37-x20) * (x37-x20) + (y37-y20) * (y37-y20) - 1.0;
subject to d37_21:
	0 >= (x37-x21) * (x37-x21) + (y37-y21) * (y37-y21) - 1.0;
subject to d37_22:
	0 >= (x37-x22) * (x37-x22) + (y37-y22) * (y37-y22) - 1.0;
subject to d37_23:
	0 >= (x37-x23) * (x37-x23) + (y37-y23) * (y37-y23) - 1.0;
subject to d37_24:
	0 >= (x37-x24) * (x37-x24) + (y37-y24) * (y37-y24) - 1.0;
subject to d37_25:
	0 >= (x37-x25) * (x37-x25) + (y37-y25) * (y37-y25) - 1.0;
subject to d37_26:
	0 >= (x37-x26) * (x37-x26) + (y37-y26) * (y37-y26) - 1.0;
subject to d37_27:
	0 >= (x37-x27) * (x37-x27) + (y37-y27) * (y37-y27) - 1.0;
subject to d37_28:
	0 >= (x37-x28) * (x37-x28) + (y37-y28) * (y37-y28) - 1.0;
subject to d37_29:
	0 >= (x37-x29) * (x37-x29) + (y37-y29) * (y37-y29) - 1.0;
subject to d37_30:
	0 >= (x37-x30) * (x37-x30) + (y37-y30) * (y37-y30) - 1.0;
subject to d37_31:
	0 >= (x37-x31) * (x37-x31) + (y37-y31) * (y37-y31) - 1.0;
subject to d37_32:
	0 >= (x37-x32) * (x37-x32) + (y37-y32) * (y37-y32) - 1.0;
subject to d37_33:
	0 >= (x37-x33) * (x37-x33) + (y37-y33) * (y37-y33) - 1.0;
subject to d37_34:
	0 >= (x37-x34) * (x37-x34) + (y37-y34) * (y37-y34) - 1.0;
subject to d37_35:
	0 >= (x37-x35) * (x37-x35) + (y37-y35) * (y37-y35) - 1.0;
subject to d37_36:
	0 >= (x37-x36) * (x37-x36) + (y37-y36) * (y37-y36) - 1.0;
subject to d38_0:
	0 >= (x38-x0) * (x38-x0) + (y38-y0) * (y38-y0) - 1.0;
subject to d38_1:
	0 >= (x38-x1) * (x38-x1) + (y38-y1) * (y38-y1) - 1.0;
subject to d38_2:
	0 >= (x38-x2) * (x38-x2) + (y38-y2) * (y38-y2) - 1.0;
subject to d38_3:
	0 >= (x38-x3) * (x38-x3) + (y38-y3) * (y38-y3) - 1.0;
subject to d38_4:
	0 >= (x38-x4) * (x38-x4) + (y38-y4) * (y38-y4) - 1.0;
subject to d38_5:
	0 >= (x38-x5) * (x38-x5) + (y38-y5) * (y38-y5) - 1.0;
subject to d38_6:
	0 >= (x38-x6) * (x38-x6) + (y38-y6) * (y38-y6) - 1.0;
subject to d38_7:
	0 >= (x38-x7) * (x38-x7) + (y38-y7) * (y38-y7) - 1.0;
subject to d38_8:
	0 >= (x38-x8) * (x38-x8) + (y38-y8) * (y38-y8) - 1.0;
subject to d38_9:
	0 >= (x38-x9) * (x38-x9) + (y38-y9) * (y38-y9) - 1.0;
subject to d38_10:
	0 >= (x38-x10) * (x38-x10) + (y38-y10) * (y38-y10) - 1.0;
subject to d38_11:
	0 >= (x38-x11) * (x38-x11) + (y38-y11) * (y38-y11) - 1.0;
subject to d38_12:
	0 >= (x38-x12) * (x38-x12) + (y38-y12) * (y38-y12) - 1.0;
subject to d38_13:
	0 >= (x38-x13) * (x38-x13) + (y38-y13) * (y38-y13) - 1.0;
subject to d38_14:
	0 >= (x38-x14) * (x38-x14) + (y38-y14) * (y38-y14) - 1.0;
subject to d38_15:
	0 >= (x38-x15) * (x38-x15) + (y38-y15) * (y38-y15) - 1.0;
subject to d38_16:
	0 >= (x38-x16) * (x38-x16) + (y38-y16) * (y38-y16) - 1.0;
subject to d38_17:
	0 >= (x38-x17) * (x38-x17) + (y38-y17) * (y38-y17) - 1.0;
subject to d38_18:
	0 >= (x38-x18) * (x38-x18) + (y38-y18) * (y38-y18) - 1.0;
subject to d38_19:
	0 >= (x38-x19) * (x38-x19) + (y38-y19) * (y38-y19) - 1.0;
subject to d38_20:
	0 >= (x38-x20) * (x38-x20) + (y38-y20) * (y38-y20) - 1.0;
subject to d38_21:
	0 >= (x38-x21) * (x38-x21) + (y38-y21) * (y38-y21) - 1.0;
subject to d38_22:
	0 >= (x38-x22) * (x38-x22) + (y38-y22) * (y38-y22) - 1.0;
subject to d38_23:
	0 >= (x38-x23) * (x38-x23) + (y38-y23) * (y38-y23) - 1.0;
subject to d38_24:
	0 >= (x38-x24) * (x38-x24) + (y38-y24) * (y38-y24) - 1.0;
subject to d38_25:
	0 >= (x38-x25) * (x38-x25) + (y38-y25) * (y38-y25) - 1.0;
subject to d38_26:
	0 >= (x38-x26) * (x38-x26) + (y38-y26) * (y38-y26) - 1.0;
subject to d38_27:
	0 >= (x38-x27) * (x38-x27) + (y38-y27) * (y38-y27) - 1.0;
subject to d38_28:
	0 >= (x38-x28) * (x38-x28) + (y38-y28) * (y38-y28) - 1.0;
subject to d38_29:
	0 >= (x38-x29) * (x38-x29) + (y38-y29) * (y38-y29) - 1.0;
subject to d38_30:
	0 >= (x38-x30) * (x38-x30) + (y38-y30) * (y38-y30) - 1.0;
subject to d38_31:
	0 >= (x38-x31) * (x38-x31) + (y38-y31) * (y38-y31) - 1.0;
subject to d38_32:
	0 >= (x38-x32) * (x38-x32) + (y38-y32) * (y38-y32) - 1.0;
subject to d38_33:
	0 >= (x38-x33) * (x38-x33) + (y38-y33) * (y38-y33) - 1.0;
subject to d38_34:
	0 >= (x38-x34) * (x38-x34) + (y38-y34) * (y38-y34) - 1.0;
subject to d38_35:
	0 >= (x38-x35) * (x38-x35) + (y38-y35) * (y38-y35) - 1.0;
subject to d38_36:
	0 >= (x38-x36) * (x38-x36) + (y38-y36) * (y38-y36) - 1.0;
subject to d38_37:
	0 >= (x38-x37) * (x38-x37) + (y38-y37) * (y38-y37) - 1.0;
subject to d39_0:
	0 >= (x39-x0) * (x39-x0) + (y39-y0) * (y39-y0) - 1.0;
subject to d39_1:
	0 >= (x39-x1) * (x39-x1) + (y39-y1) * (y39-y1) - 1.0;
subject to d39_2:
	0 >= (x39-x2) * (x39-x2) + (y39-y2) * (y39-y2) - 1.0;
subject to d39_3:
	0 >= (x39-x3) * (x39-x3) + (y39-y3) * (y39-y3) - 1.0;
subject to d39_4:
	0 >= (x39-x4) * (x39-x4) + (y39-y4) * (y39-y4) - 1.0;
subject to d39_5:
	0 >= (x39-x5) * (x39-x5) + (y39-y5) * (y39-y5) - 1.0;
subject to d39_6:
	0 >= (x39-x6) * (x39-x6) + (y39-y6) * (y39-y6) - 1.0;
subject to d39_7:
	0 >= (x39-x7) * (x39-x7) + (y39-y7) * (y39-y7) - 1.0;
subject to d39_8:
	0 >= (x39-x8) * (x39-x8) + (y39-y8) * (y39-y8) - 1.0;
subject to d39_9:
	0 >= (x39-x9) * (x39-x9) + (y39-y9) * (y39-y9) - 1.0;
subject to d39_10:
	0 >= (x39-x10) * (x39-x10) + (y39-y10) * (y39-y10) - 1.0;
subject to d39_11:
	0 >= (x39-x11) * (x39-x11) + (y39-y11) * (y39-y11) - 1.0;
subject to d39_12:
	0 >= (x39-x12) * (x39-x12) + (y39-y12) * (y39-y12) - 1.0;
subject to d39_13:
	0 >= (x39-x13) * (x39-x13) + (y39-y13) * (y39-y13) - 1.0;
subject to d39_14:
	0 >= (x39-x14) * (x39-x14) + (y39-y14) * (y39-y14) - 1.0;
subject to d39_15:
	0 >= (x39-x15) * (x39-x15) + (y39-y15) * (y39-y15) - 1.0;
subject to d39_16:
	0 >= (x39-x16) * (x39-x16) + (y39-y16) * (y39-y16) - 1.0;
subject to d39_17:
	0 >= (x39-x17) * (x39-x17) + (y39-y17) * (y39-y17) - 1.0;
subject to d39_18:
	0 >= (x39-x18) * (x39-x18) + (y39-y18) * (y39-y18) - 1.0;
subject to d39_19:
	0 >= (x39-x19) * (x39-x19) + (y39-y19) * (y39-y19) - 1.0;
subject to d39_20:
	0 >= (x39-x20) * (x39-x20) + (y39-y20) * (y39-y20) - 1.0;
subject to d39_21:
	0 >= (x39-x21) * (x39-x21) + (y39-y21) * (y39-y21) - 1.0;
subject to d39_22:
	0 >= (x39-x22) * (x39-x22) + (y39-y22) * (y39-y22) - 1.0;
subject to d39_23:
	0 >= (x39-x23) * (x39-x23) + (y39-y23) * (y39-y23) - 1.0;
subject to d39_24:
	0 >= (x39-x24) * (x39-x24) + (y39-y24) * (y39-y24) - 1.0;
subject to d39_25:
	0 >= (x39-x25) * (x39-x25) + (y39-y25) * (y39-y25) - 1.0;
subject to d39_26:
	0 >= (x39-x26) * (x39-x26) + (y39-y26) * (y39-y26) - 1.0;
subject to d39_27:
	0 >= (x39-x27) * (x39-x27) + (y39-y27) * (y39-y27) - 1.0;
subject to d39_28:
	0 >= (x39-x28) * (x39-x28) + (y39-y28) * (y39-y28) - 1.0;
subject to d39_29:
	0 >= (x39-x29) * (x39-x29) + (y39-y29) * (y39-y29) - 1.0;
subject to d39_30:
	0 >= (x39-x30) * (x39-x30) + (y39-y30) * (y39-y30) - 1.0;
subject to d39_31:
	0 >= (x39-x31) * (x39-x31) + (y39-y31) * (y39-y31) - 1.0;
subject to d39_32:
	0 >= (x39-x32) * (x39-x32) + (y39-y32) * (y39-y32) - 1.0;
subject to d39_33:
	0 >= (x39-x33) * (x39-x33) + (y39-y33) * (y39-y33) - 1.0;
subject to d39_34:
	0 >= (x39-x34) * (x39-x34) + (y39-y34) * (y39-y34) - 1.0;
subject to d39_35:
	0 >= (x39-x35) * (x39-x35) + (y39-y35) * (y39-y35) - 1.0;
subject to d39_36:
	0 >= (x39-x36) * (x39-x36) + (y39-y36) * (y39-y36) - 1.0;
subject to d39_37:
	0 >= (x39-x37) * (x39-x37) + (y39-y37) * (y39-y37) - 1.0;
subject to d39_38:
	0 >= (x39-x38) * (x39-x38) + (y39-y38) * (y39-y38) - 1.0;
subject to d40_0:
	0 >= (x40-x0) * (x40-x0) + (y40-y0) * (y40-y0) - 1.0;
subject to d40_1:
	0 >= (x40-x1) * (x40-x1) + (y40-y1) * (y40-y1) - 1.0;
subject to d40_2:
	0 >= (x40-x2) * (x40-x2) + (y40-y2) * (y40-y2) - 1.0;
subject to d40_3:
	0 >= (x40-x3) * (x40-x3) + (y40-y3) * (y40-y3) - 1.0;
subject to d40_4:
	0 >= (x40-x4) * (x40-x4) + (y40-y4) * (y40-y4) - 1.0;
subject to d40_5:
	0 >= (x40-x5) * (x40-x5) + (y40-y5) * (y40-y5) - 1.0;
subject to d40_6:
	0 >= (x40-x6) * (x40-x6) + (y40-y6) * (y40-y6) - 1.0;
subject to d40_7:
	0 >= (x40-x7) * (x40-x7) + (y40-y7) * (y40-y7) - 1.0;
subject to d40_8:
	0 >= (x40-x8) * (x40-x8) + (y40-y8) * (y40-y8) - 1.0;
subject to d40_9:
	0 >= (x40-x9) * (x40-x9) + (y40-y9) * (y40-y9) - 1.0;
subject to d40_10:
	0 >= (x40-x10) * (x40-x10) + (y40-y10) * (y40-y10) - 1.0;
subject to d40_11:
	0 >= (x40-x11) * (x40-x11) + (y40-y11) * (y40-y11) - 1.0;
subject to d40_12:
	0 >= (x40-x12) * (x40-x12) + (y40-y12) * (y40-y12) - 1.0;
subject to d40_13:
	0 >= (x40-x13) * (x40-x13) + (y40-y13) * (y40-y13) - 1.0;
subject to d40_14:
	0 >= (x40-x14) * (x40-x14) + (y40-y14) * (y40-y14) - 1.0;
subject to d40_15:
	0 >= (x40-x15) * (x40-x15) + (y40-y15) * (y40-y15) - 1.0;
subject to d40_16:
	0 >= (x40-x16) * (x40-x16) + (y40-y16) * (y40-y16) - 1.0;
subject to d40_17:
	0 >= (x40-x17) * (x40-x17) + (y40-y17) * (y40-y17) - 1.0;
subject to d40_18:
	0 >= (x40-x18) * (x40-x18) + (y40-y18) * (y40-y18) - 1.0;
subject to d40_19:
	0 >= (x40-x19) * (x40-x19) + (y40-y19) * (y40-y19) - 1.0;
subject to d40_20:
	0 >= (x40-x20) * (x40-x20) + (y40-y20) * (y40-y20) - 1.0;
subject to d40_21:
	0 >= (x40-x21) * (x40-x21) + (y40-y21) * (y40-y21) - 1.0;
subject to d40_22:
	0 >= (x40-x22) * (x40-x22) + (y40-y22) * (y40-y22) - 1.0;
subject to d40_23:
	0 >= (x40-x23) * (x40-x23) + (y40-y23) * (y40-y23) - 1.0;
subject to d40_24:
	0 >= (x40-x24) * (x40-x24) + (y40-y24) * (y40-y24) - 1.0;
subject to d40_25:
	0 >= (x40-x25) * (x40-x25) + (y40-y25) * (y40-y25) - 1.0;
subject to d40_26:
	0 >= (x40-x26) * (x40-x26) + (y40-y26) * (y40-y26) - 1.0;
subject to d40_27:
	0 >= (x40-x27) * (x40-x27) + (y40-y27) * (y40-y27) - 1.0;
subject to d40_28:
	0 >= (x40-x28) * (x40-x28) + (y40-y28) * (y40-y28) - 1.0;
subject to d40_29:
	0 >= (x40-x29) * (x40-x29) + (y40-y29) * (y40-y29) - 1.0;
subject to d40_30:
	0 >= (x40-x30) * (x40-x30) + (y40-y30) * (y40-y30) - 1.0;
subject to d40_31:
	0 >= (x40-x31) * (x40-x31) + (y40-y31) * (y40-y31) - 1.0;
subject to d40_32:
	0 >= (x40-x32) * (x40-x32) + (y40-y32) * (y40-y32) - 1.0;
subject to d40_33:
	0 >= (x40-x33) * (x40-x33) + (y40-y33) * (y40-y33) - 1.0;
subject to d40_34:
	0 >= (x40-x34) * (x40-x34) + (y40-y34) * (y40-y34) - 1.0;
subject to d40_35:
	0 >= (x40-x35) * (x40-x35) + (y40-y35) * (y40-y35) - 1.0;
subject to d40_36:
	0 >= (x40-x36) * (x40-x36) + (y40-y36) * (y40-y36) - 1.0;
subject to d40_37:
	0 >= (x40-x37) * (x40-x37) + (y40-y37) * (y40-y37) - 1.0;
subject to d40_38:
	0 >= (x40-x38) * (x40-x38) + (y40-y38) * (y40-y38) - 1.0;
subject to d40_39:
	0 >= (x40-x39) * (x40-x39) + (y40-y39) * (y40-y39) - 1.0;
subject to d41_0:
	0 >= (x41-x0) * (x41-x0) + (y41-y0) * (y41-y0) - 1.0;
subject to d41_1:
	0 >= (x41-x1) * (x41-x1) + (y41-y1) * (y41-y1) - 1.0;
subject to d41_2:
	0 >= (x41-x2) * (x41-x2) + (y41-y2) * (y41-y2) - 1.0;
subject to d41_3:
	0 >= (x41-x3) * (x41-x3) + (y41-y3) * (y41-y3) - 1.0;
subject to d41_4:
	0 >= (x41-x4) * (x41-x4) + (y41-y4) * (y41-y4) - 1.0;
subject to d41_5:
	0 >= (x41-x5) * (x41-x5) + (y41-y5) * (y41-y5) - 1.0;
subject to d41_6:
	0 >= (x41-x6) * (x41-x6) + (y41-y6) * (y41-y6) - 1.0;
subject to d41_7:
	0 >= (x41-x7) * (x41-x7) + (y41-y7) * (y41-y7) - 1.0;
subject to d41_8:
	0 >= (x41-x8) * (x41-x8) + (y41-y8) * (y41-y8) - 1.0;
subject to d41_9:
	0 >= (x41-x9) * (x41-x9) + (y41-y9) * (y41-y9) - 1.0;
subject to d41_10:
	0 >= (x41-x10) * (x41-x10) + (y41-y10) * (y41-y10) - 1.0;
subject to d41_11:
	0 >= (x41-x11) * (x41-x11) + (y41-y11) * (y41-y11) - 1.0;
subject to d41_12:
	0 >= (x41-x12) * (x41-x12) + (y41-y12) * (y41-y12) - 1.0;
subject to d41_13:
	0 >= (x41-x13) * (x41-x13) + (y41-y13) * (y41-y13) - 1.0;
subject to d41_14:
	0 >= (x41-x14) * (x41-x14) + (y41-y14) * (y41-y14) - 1.0;
subject to d41_15:
	0 >= (x41-x15) * (x41-x15) + (y41-y15) * (y41-y15) - 1.0;
subject to d41_16:
	0 >= (x41-x16) * (x41-x16) + (y41-y16) * (y41-y16) - 1.0;
subject to d41_17:
	0 >= (x41-x17) * (x41-x17) + (y41-y17) * (y41-y17) - 1.0;
subject to d41_18:
	0 >= (x41-x18) * (x41-x18) + (y41-y18) * (y41-y18) - 1.0;
subject to d41_19:
	0 >= (x41-x19) * (x41-x19) + (y41-y19) * (y41-y19) - 1.0;
subject to d41_20:
	0 >= (x41-x20) * (x41-x20) + (y41-y20) * (y41-y20) - 1.0;
subject to d41_21:
	0 >= (x41-x21) * (x41-x21) + (y41-y21) * (y41-y21) - 1.0;
subject to d41_22:
	0 >= (x41-x22) * (x41-x22) + (y41-y22) * (y41-y22) - 1.0;
subject to d41_23:
	0 >= (x41-x23) * (x41-x23) + (y41-y23) * (y41-y23) - 1.0;
subject to d41_24:
	0 >= (x41-x24) * (x41-x24) + (y41-y24) * (y41-y24) - 1.0;
subject to d41_25:
	0 >= (x41-x25) * (x41-x25) + (y41-y25) * (y41-y25) - 1.0;
subject to d41_26:
	0 >= (x41-x26) * (x41-x26) + (y41-y26) * (y41-y26) - 1.0;
subject to d41_27:
	0 >= (x41-x27) * (x41-x27) + (y41-y27) * (y41-y27) - 1.0;
subject to d41_28:
	0 >= (x41-x28) * (x41-x28) + (y41-y28) * (y41-y28) - 1.0;
subject to d41_29:
	0 >= (x41-x29) * (x41-x29) + (y41-y29) * (y41-y29) - 1.0;
subject to d41_30:
	0 >= (x41-x30) * (x41-x30) + (y41-y30) * (y41-y30) - 1.0;
subject to d41_31:
	0 >= (x41-x31) * (x41-x31) + (y41-y31) * (y41-y31) - 1.0;
subject to d41_32:
	0 >= (x41-x32) * (x41-x32) + (y41-y32) * (y41-y32) - 1.0;
subject to d41_33:
	0 >= (x41-x33) * (x41-x33) + (y41-y33) * (y41-y33) - 1.0;
subject to d41_34:
	0 >= (x41-x34) * (x41-x34) + (y41-y34) * (y41-y34) - 1.0;
subject to d41_35:
	0 >= (x41-x35) * (x41-x35) + (y41-y35) * (y41-y35) - 1.0;
subject to d41_36:
	0 >= (x41-x36) * (x41-x36) + (y41-y36) * (y41-y36) - 1.0;
subject to d41_37:
	0 >= (x41-x37) * (x41-x37) + (y41-y37) * (y41-y37) - 1.0;
subject to d41_38:
	0 >= (x41-x38) * (x41-x38) + (y41-y38) * (y41-y38) - 1.0;
subject to d41_39:
	0 >= (x41-x39) * (x41-x39) + (y41-y39) * (y41-y39) - 1.0;
subject to d41_40:
	0 >= (x41-x40) * (x41-x40) + (y41-y40) * (y41-y40) - 1.0;
subject to d42_0:
	0 >= (x42-x0) * (x42-x0) + (y42-y0) * (y42-y0) - 1.0;
subject to d42_1:
	0 >= (x42-x1) * (x42-x1) + (y42-y1) * (y42-y1) - 1.0;
subject to d42_2:
	0 >= (x42-x2) * (x42-x2) + (y42-y2) * (y42-y2) - 1.0;
subject to d42_3:
	0 >= (x42-x3) * (x42-x3) + (y42-y3) * (y42-y3) - 1.0;
subject to d42_4:
	0 >= (x42-x4) * (x42-x4) + (y42-y4) * (y42-y4) - 1.0;
subject to d42_5:
	0 >= (x42-x5) * (x42-x5) + (y42-y5) * (y42-y5) - 1.0;
subject to d42_6:
	0 >= (x42-x6) * (x42-x6) + (y42-y6) * (y42-y6) - 1.0;
subject to d42_7:
	0 >= (x42-x7) * (x42-x7) + (y42-y7) * (y42-y7) - 1.0;
subject to d42_8:
	0 >= (x42-x8) * (x42-x8) + (y42-y8) * (y42-y8) - 1.0;
subject to d42_9:
	0 >= (x42-x9) * (x42-x9) + (y42-y9) * (y42-y9) - 1.0;
subject to d42_10:
	0 >= (x42-x10) * (x42-x10) + (y42-y10) * (y42-y10) - 1.0;
subject to d42_11:
	0 >= (x42-x11) * (x42-x11) + (y42-y11) * (y42-y11) - 1.0;
subject to d42_12:
	0 >= (x42-x12) * (x42-x12) + (y42-y12) * (y42-y12) - 1.0;
subject to d42_13:
	0 >= (x42-x13) * (x42-x13) + (y42-y13) * (y42-y13) - 1.0;
subject to d42_14:
	0 >= (x42-x14) * (x42-x14) + (y42-y14) * (y42-y14) - 1.0;
subject to d42_15:
	0 >= (x42-x15) * (x42-x15) + (y42-y15) * (y42-y15) - 1.0;
subject to d42_16:
	0 >= (x42-x16) * (x42-x16) + (y42-y16) * (y42-y16) - 1.0;
subject to d42_17:
	0 >= (x42-x17) * (x42-x17) + (y42-y17) * (y42-y17) - 1.0;
subject to d42_18:
	0 >= (x42-x18) * (x42-x18) + (y42-y18) * (y42-y18) - 1.0;
subject to d42_19:
	0 >= (x42-x19) * (x42-x19) + (y42-y19) * (y42-y19) - 1.0;
subject to d42_20:
	0 >= (x42-x20) * (x42-x20) + (y42-y20) * (y42-y20) - 1.0;
subject to d42_21:
	0 >= (x42-x21) * (x42-x21) + (y42-y21) * (y42-y21) - 1.0;
subject to d42_22:
	0 >= (x42-x22) * (x42-x22) + (y42-y22) * (y42-y22) - 1.0;
subject to d42_23:
	0 >= (x42-x23) * (x42-x23) + (y42-y23) * (y42-y23) - 1.0;
subject to d42_24:
	0 >= (x42-x24) * (x42-x24) + (y42-y24) * (y42-y24) - 1.0;
subject to d42_25:
	0 >= (x42-x25) * (x42-x25) + (y42-y25) * (y42-y25) - 1.0;
subject to d42_26:
	0 >= (x42-x26) * (x42-x26) + (y42-y26) * (y42-y26) - 1.0;
subject to d42_27:
	0 >= (x42-x27) * (x42-x27) + (y42-y27) * (y42-y27) - 1.0;
subject to d42_28:
	0 >= (x42-x28) * (x42-x28) + (y42-y28) * (y42-y28) - 1.0;
subject to d42_29:
	0 >= (x42-x29) * (x42-x29) + (y42-y29) * (y42-y29) - 1.0;
subject to d42_30:
	0 >= (x42-x30) * (x42-x30) + (y42-y30) * (y42-y30) - 1.0;
subject to d42_31:
	0 >= (x42-x31) * (x42-x31) + (y42-y31) * (y42-y31) - 1.0;
subject to d42_32:
	0 >= (x42-x32) * (x42-x32) + (y42-y32) * (y42-y32) - 1.0;
subject to d42_33:
	0 >= (x42-x33) * (x42-x33) + (y42-y33) * (y42-y33) - 1.0;
subject to d42_34:
	0 >= (x42-x34) * (x42-x34) + (y42-y34) * (y42-y34) - 1.0;
subject to d42_35:
	0 >= (x42-x35) * (x42-x35) + (y42-y35) * (y42-y35) - 1.0;
subject to d42_36:
	0 >= (x42-x36) * (x42-x36) + (y42-y36) * (y42-y36) - 1.0;
subject to d42_37:
	0 >= (x42-x37) * (x42-x37) + (y42-y37) * (y42-y37) - 1.0;
subject to d42_38:
	0 >= (x42-x38) * (x42-x38) + (y42-y38) * (y42-y38) - 1.0;
subject to d42_39:
	0 >= (x42-x39) * (x42-x39) + (y42-y39) * (y42-y39) - 1.0;
subject to d42_40:
	0 >= (x42-x40) * (x42-x40) + (y42-y40) * (y42-y40) - 1.0;
subject to d42_41:
	0 >= (x42-x41) * (x42-x41) + (y42-y41) * (y42-y41) - 1.0;
subject to d43_0:
	0 >= (x43-x0) * (x43-x0) + (y43-y0) * (y43-y0) - 1.0;
subject to d43_1:
	0 >= (x43-x1) * (x43-x1) + (y43-y1) * (y43-y1) - 1.0;
subject to d43_2:
	0 >= (x43-x2) * (x43-x2) + (y43-y2) * (y43-y2) - 1.0;
subject to d43_3:
	0 >= (x43-x3) * (x43-x3) + (y43-y3) * (y43-y3) - 1.0;
subject to d43_4:
	0 >= (x43-x4) * (x43-x4) + (y43-y4) * (y43-y4) - 1.0;
subject to d43_5:
	0 >= (x43-x5) * (x43-x5) + (y43-y5) * (y43-y5) - 1.0;
subject to d43_6:
	0 >= (x43-x6) * (x43-x6) + (y43-y6) * (y43-y6) - 1.0;
subject to d43_7:
	0 >= (x43-x7) * (x43-x7) + (y43-y7) * (y43-y7) - 1.0;
subject to d43_8:
	0 >= (x43-x8) * (x43-x8) + (y43-y8) * (y43-y8) - 1.0;
subject to d43_9:
	0 >= (x43-x9) * (x43-x9) + (y43-y9) * (y43-y9) - 1.0;
subject to d43_10:
	0 >= (x43-x10) * (x43-x10) + (y43-y10) * (y43-y10) - 1.0;
subject to d43_11:
	0 >= (x43-x11) * (x43-x11) + (y43-y11) * (y43-y11) - 1.0;
subject to d43_12:
	0 >= (x43-x12) * (x43-x12) + (y43-y12) * (y43-y12) - 1.0;
subject to d43_13:
	0 >= (x43-x13) * (x43-x13) + (y43-y13) * (y43-y13) - 1.0;
subject to d43_14:
	0 >= (x43-x14) * (x43-x14) + (y43-y14) * (y43-y14) - 1.0;
subject to d43_15:
	0 >= (x43-x15) * (x43-x15) + (y43-y15) * (y43-y15) - 1.0;
subject to d43_16:
	0 >= (x43-x16) * (x43-x16) + (y43-y16) * (y43-y16) - 1.0;
subject to d43_17:
	0 >= (x43-x17) * (x43-x17) + (y43-y17) * (y43-y17) - 1.0;
subject to d43_18:
	0 >= (x43-x18) * (x43-x18) + (y43-y18) * (y43-y18) - 1.0;
subject to d43_19:
	0 >= (x43-x19) * (x43-x19) + (y43-y19) * (y43-y19) - 1.0;
subject to d43_20:
	0 >= (x43-x20) * (x43-x20) + (y43-y20) * (y43-y20) - 1.0;
subject to d43_21:
	0 >= (x43-x21) * (x43-x21) + (y43-y21) * (y43-y21) - 1.0;
subject to d43_22:
	0 >= (x43-x22) * (x43-x22) + (y43-y22) * (y43-y22) - 1.0;
subject to d43_23:
	0 >= (x43-x23) * (x43-x23) + (y43-y23) * (y43-y23) - 1.0;
subject to d43_24:
	0 >= (x43-x24) * (x43-x24) + (y43-y24) * (y43-y24) - 1.0;
subject to d43_25:
	0 >= (x43-x25) * (x43-x25) + (y43-y25) * (y43-y25) - 1.0;
subject to d43_26:
	0 >= (x43-x26) * (x43-x26) + (y43-y26) * (y43-y26) - 1.0;
subject to d43_27:
	0 >= (x43-x27) * (x43-x27) + (y43-y27) * (y43-y27) - 1.0;
subject to d43_28:
	0 >= (x43-x28) * (x43-x28) + (y43-y28) * (y43-y28) - 1.0;
subject to d43_29:
	0 >= (x43-x29) * (x43-x29) + (y43-y29) * (y43-y29) - 1.0;
subject to d43_30:
	0 >= (x43-x30) * (x43-x30) + (y43-y30) * (y43-y30) - 1.0;
subject to d43_31:
	0 >= (x43-x31) * (x43-x31) + (y43-y31) * (y43-y31) - 1.0;
subject to d43_32:
	0 >= (x43-x32) * (x43-x32) + (y43-y32) * (y43-y32) - 1.0;
subject to d43_33:
	0 >= (x43-x33) * (x43-x33) + (y43-y33) * (y43-y33) - 1.0;
subject to d43_34:
	0 >= (x43-x34) * (x43-x34) + (y43-y34) * (y43-y34) - 1.0;
subject to d43_35:
	0 >= (x43-x35) * (x43-x35) + (y43-y35) * (y43-y35) - 1.0;
subject to d43_36:
	0 >= (x43-x36) * (x43-x36) + (y43-y36) * (y43-y36) - 1.0;
subject to d43_37:
	0 >= (x43-x37) * (x43-x37) + (y43-y37) * (y43-y37) - 1.0;
subject to d43_38:
	0 >= (x43-x38) * (x43-x38) + (y43-y38) * (y43-y38) - 1.0;
subject to d43_39:
	0 >= (x43-x39) * (x43-x39) + (y43-y39) * (y43-y39) - 1.0;
subject to d43_40:
	0 >= (x43-x40) * (x43-x40) + (y43-y40) * (y43-y40) - 1.0;
subject to d43_41:
	0 >= (x43-x41) * (x43-x41) + (y43-y41) * (y43-y41) - 1.0;
subject to d43_42:
	0 >= (x43-x42) * (x43-x42) + (y43-y42) * (y43-y42) - 1.0;
subject to d44_0:
	0 >= (x44-x0) * (x44-x0) + (y44-y0) * (y44-y0) - 1.0;
subject to d44_1:
	0 >= (x44-x1) * (x44-x1) + (y44-y1) * (y44-y1) - 1.0;
subject to d44_2:
	0 >= (x44-x2) * (x44-x2) + (y44-y2) * (y44-y2) - 1.0;
subject to d44_3:
	0 >= (x44-x3) * (x44-x3) + (y44-y3) * (y44-y3) - 1.0;
subject to d44_4:
	0 >= (x44-x4) * (x44-x4) + (y44-y4) * (y44-y4) - 1.0;
subject to d44_5:
	0 >= (x44-x5) * (x44-x5) + (y44-y5) * (y44-y5) - 1.0;
subject to d44_6:
	0 >= (x44-x6) * (x44-x6) + (y44-y6) * (y44-y6) - 1.0;
subject to d44_7:
	0 >= (x44-x7) * (x44-x7) + (y44-y7) * (y44-y7) - 1.0;
subject to d44_8:
	0 >= (x44-x8) * (x44-x8) + (y44-y8) * (y44-y8) - 1.0;
subject to d44_9:
	0 >= (x44-x9) * (x44-x9) + (y44-y9) * (y44-y9) - 1.0;
subject to d44_10:
	0 >= (x44-x10) * (x44-x10) + (y44-y10) * (y44-y10) - 1.0;
subject to d44_11:
	0 >= (x44-x11) * (x44-x11) + (y44-y11) * (y44-y11) - 1.0;
subject to d44_12:
	0 >= (x44-x12) * (x44-x12) + (y44-y12) * (y44-y12) - 1.0;
subject to d44_13:
	0 >= (x44-x13) * (x44-x13) + (y44-y13) * (y44-y13) - 1.0;
subject to d44_14:
	0 >= (x44-x14) * (x44-x14) + (y44-y14) * (y44-y14) - 1.0;
subject to d44_15:
	0 >= (x44-x15) * (x44-x15) + (y44-y15) * (y44-y15) - 1.0;
subject to d44_16:
	0 >= (x44-x16) * (x44-x16) + (y44-y16) * (y44-y16) - 1.0;
subject to d44_17:
	0 >= (x44-x17) * (x44-x17) + (y44-y17) * (y44-y17) - 1.0;
subject to d44_18:
	0 >= (x44-x18) * (x44-x18) + (y44-y18) * (y44-y18) - 1.0;
subject to d44_19:
	0 >= (x44-x19) * (x44-x19) + (y44-y19) * (y44-y19) - 1.0;
subject to d44_20:
	0 >= (x44-x20) * (x44-x20) + (y44-y20) * (y44-y20) - 1.0;
subject to d44_21:
	0 >= (x44-x21) * (x44-x21) + (y44-y21) * (y44-y21) - 1.0;
subject to d44_22:
	0 >= (x44-x22) * (x44-x22) + (y44-y22) * (y44-y22) - 1.0;
subject to d44_23:
	0 >= (x44-x23) * (x44-x23) + (y44-y23) * (y44-y23) - 1.0;
subject to d44_24:
	0 >= (x44-x24) * (x44-x24) + (y44-y24) * (y44-y24) - 1.0;
subject to d44_25:
	0 >= (x44-x25) * (x44-x25) + (y44-y25) * (y44-y25) - 1.0;
subject to d44_26:
	0 >= (x44-x26) * (x44-x26) + (y44-y26) * (y44-y26) - 1.0;
subject to d44_27:
	0 >= (x44-x27) * (x44-x27) + (y44-y27) * (y44-y27) - 1.0;
subject to d44_28:
	0 >= (x44-x28) * (x44-x28) + (y44-y28) * (y44-y28) - 1.0;
subject to d44_29:
	0 >= (x44-x29) * (x44-x29) + (y44-y29) * (y44-y29) - 1.0;
subject to d44_30:
	0 >= (x44-x30) * (x44-x30) + (y44-y30) * (y44-y30) - 1.0;
subject to d44_31:
	0 >= (x44-x31) * (x44-x31) + (y44-y31) * (y44-y31) - 1.0;
subject to d44_32:
	0 >= (x44-x32) * (x44-x32) + (y44-y32) * (y44-y32) - 1.0;
subject to d44_33:
	0 >= (x44-x33) * (x44-x33) + (y44-y33) * (y44-y33) - 1.0;
subject to d44_34:
	0 >= (x44-x34) * (x44-x34) + (y44-y34) * (y44-y34) - 1.0;
subject to d44_35:
	0 >= (x44-x35) * (x44-x35) + (y44-y35) * (y44-y35) - 1.0;
subject to d44_36:
	0 >= (x44-x36) * (x44-x36) + (y44-y36) * (y44-y36) - 1.0;
subject to d44_37:
	0 >= (x44-x37) * (x44-x37) + (y44-y37) * (y44-y37) - 1.0;
subject to d44_38:
	0 >= (x44-x38) * (x44-x38) + (y44-y38) * (y44-y38) - 1.0;
subject to d44_39:
	0 >= (x44-x39) * (x44-x39) + (y44-y39) * (y44-y39) - 1.0;
subject to d44_40:
	0 >= (x44-x40) * (x44-x40) + (y44-y40) * (y44-y40) - 1.0;
subject to d44_41:
	0 >= (x44-x41) * (x44-x41) + (y44-y41) * (y44-y41) - 1.0;
subject to d44_42:
	0 >= (x44-x42) * (x44-x42) + (y44-y42) * (y44-y42) - 1.0;
subject to d44_43:
	0 >= (x44-x43) * (x44-x43) + (y44-y43) * (y44-y43) - 1.0;
subject to d45_0:
	0 >= (x45-x0) * (x45-x0) + (y45-y0) * (y45-y0) - 1.0;
subject to d45_1:
	0 >= (x45-x1) * (x45-x1) + (y45-y1) * (y45-y1) - 1.0;
subject to d45_2:
	0 >= (x45-x2) * (x45-x2) + (y45-y2) * (y45-y2) - 1.0;
subject to d45_3:
	0 >= (x45-x3) * (x45-x3) + (y45-y3) * (y45-y3) - 1.0;
subject to d45_4:
	0 >= (x45-x4) * (x45-x4) + (y45-y4) * (y45-y4) - 1.0;
subject to d45_5:
	0 >= (x45-x5) * (x45-x5) + (y45-y5) * (y45-y5) - 1.0;
subject to d45_6:
	0 >= (x45-x6) * (x45-x6) + (y45-y6) * (y45-y6) - 1.0;
subject to d45_7:
	0 >= (x45-x7) * (x45-x7) + (y45-y7) * (y45-y7) - 1.0;
subject to d45_8:
	0 >= (x45-x8) * (x45-x8) + (y45-y8) * (y45-y8) - 1.0;
subject to d45_9:
	0 >= (x45-x9) * (x45-x9) + (y45-y9) * (y45-y9) - 1.0;
subject to d45_10:
	0 >= (x45-x10) * (x45-x10) + (y45-y10) * (y45-y10) - 1.0;
subject to d45_11:
	0 >= (x45-x11) * (x45-x11) + (y45-y11) * (y45-y11) - 1.0;
subject to d45_12:
	0 >= (x45-x12) * (x45-x12) + (y45-y12) * (y45-y12) - 1.0;
subject to d45_13:
	0 >= (x45-x13) * (x45-x13) + (y45-y13) * (y45-y13) - 1.0;
subject to d45_14:
	0 >= (x45-x14) * (x45-x14) + (y45-y14) * (y45-y14) - 1.0;
subject to d45_15:
	0 >= (x45-x15) * (x45-x15) + (y45-y15) * (y45-y15) - 1.0;
subject to d45_16:
	0 >= (x45-x16) * (x45-x16) + (y45-y16) * (y45-y16) - 1.0;
subject to d45_17:
	0 >= (x45-x17) * (x45-x17) + (y45-y17) * (y45-y17) - 1.0;
subject to d45_18:
	0 >= (x45-x18) * (x45-x18) + (y45-y18) * (y45-y18) - 1.0;
subject to d45_19:
	0 >= (x45-x19) * (x45-x19) + (y45-y19) * (y45-y19) - 1.0;
subject to d45_20:
	0 >= (x45-x20) * (x45-x20) + (y45-y20) * (y45-y20) - 1.0;
subject to d45_21:
	0 >= (x45-x21) * (x45-x21) + (y45-y21) * (y45-y21) - 1.0;
subject to d45_22:
	0 >= (x45-x22) * (x45-x22) + (y45-y22) * (y45-y22) - 1.0;
subject to d45_23:
	0 >= (x45-x23) * (x45-x23) + (y45-y23) * (y45-y23) - 1.0;
subject to d45_24:
	0 >= (x45-x24) * (x45-x24) + (y45-y24) * (y45-y24) - 1.0;
subject to d45_25:
	0 >= (x45-x25) * (x45-x25) + (y45-y25) * (y45-y25) - 1.0;
subject to d45_26:
	0 >= (x45-x26) * (x45-x26) + (y45-y26) * (y45-y26) - 1.0;
subject to d45_27:
	0 >= (x45-x27) * (x45-x27) + (y45-y27) * (y45-y27) - 1.0;
subject to d45_28:
	0 >= (x45-x28) * (x45-x28) + (y45-y28) * (y45-y28) - 1.0;
subject to d45_29:
	0 >= (x45-x29) * (x45-x29) + (y45-y29) * (y45-y29) - 1.0;
subject to d45_30:
	0 >= (x45-x30) * (x45-x30) + (y45-y30) * (y45-y30) - 1.0;
subject to d45_31:
	0 >= (x45-x31) * (x45-x31) + (y45-y31) * (y45-y31) - 1.0;
subject to d45_32:
	0 >= (x45-x32) * (x45-x32) + (y45-y32) * (y45-y32) - 1.0;
subject to d45_33:
	0 >= (x45-x33) * (x45-x33) + (y45-y33) * (y45-y33) - 1.0;
subject to d45_34:
	0 >= (x45-x34) * (x45-x34) + (y45-y34) * (y45-y34) - 1.0;
subject to d45_35:
	0 >= (x45-x35) * (x45-x35) + (y45-y35) * (y45-y35) - 1.0;
subject to d45_36:
	0 >= (x45-x36) * (x45-x36) + (y45-y36) * (y45-y36) - 1.0;
subject to d45_37:
	0 >= (x45-x37) * (x45-x37) + (y45-y37) * (y45-y37) - 1.0;
subject to d45_38:
	0 >= (x45-x38) * (x45-x38) + (y45-y38) * (y45-y38) - 1.0;
subject to d45_39:
	0 >= (x45-x39) * (x45-x39) + (y45-y39) * (y45-y39) - 1.0;
subject to d45_40:
	0 >= (x45-x40) * (x45-x40) + (y45-y40) * (y45-y40) - 1.0;
subject to d45_41:
	0 >= (x45-x41) * (x45-x41) + (y45-y41) * (y45-y41) - 1.0;
subject to d45_42:
	0 >= (x45-x42) * (x45-x42) + (y45-y42) * (y45-y42) - 1.0;
subject to d45_43:
	0 >= (x45-x43) * (x45-x43) + (y45-y43) * (y45-y43) - 1.0;
subject to d45_44:
	0 >= (x45-x44) * (x45-x44) + (y45-y44) * (y45-y44) - 1.0;
subject to d46_0:
	0 >= (x46-x0) * (x46-x0) + (y46-y0) * (y46-y0) - 1.0;
subject to d46_1:
	0 >= (x46-x1) * (x46-x1) + (y46-y1) * (y46-y1) - 1.0;
subject to d46_2:
	0 >= (x46-x2) * (x46-x2) + (y46-y2) * (y46-y2) - 1.0;
subject to d46_3:
	0 >= (x46-x3) * (x46-x3) + (y46-y3) * (y46-y3) - 1.0;
subject to d46_4:
	0 >= (x46-x4) * (x46-x4) + (y46-y4) * (y46-y4) - 1.0;
subject to d46_5:
	0 >= (x46-x5) * (x46-x5) + (y46-y5) * (y46-y5) - 1.0;
subject to d46_6:
	0 >= (x46-x6) * (x46-x6) + (y46-y6) * (y46-y6) - 1.0;
subject to d46_7:
	0 >= (x46-x7) * (x46-x7) + (y46-y7) * (y46-y7) - 1.0;
subject to d46_8:
	0 >= (x46-x8) * (x46-x8) + (y46-y8) * (y46-y8) - 1.0;
subject to d46_9:
	0 >= (x46-x9) * (x46-x9) + (y46-y9) * (y46-y9) - 1.0;
subject to d46_10:
	0 >= (x46-x10) * (x46-x10) + (y46-y10) * (y46-y10) - 1.0;
subject to d46_11:
	0 >= (x46-x11) * (x46-x11) + (y46-y11) * (y46-y11) - 1.0;
subject to d46_12:
	0 >= (x46-x12) * (x46-x12) + (y46-y12) * (y46-y12) - 1.0;
subject to d46_13:
	0 >= (x46-x13) * (x46-x13) + (y46-y13) * (y46-y13) - 1.0;
subject to d46_14:
	0 >= (x46-x14) * (x46-x14) + (y46-y14) * (y46-y14) - 1.0;
subject to d46_15:
	0 >= (x46-x15) * (x46-x15) + (y46-y15) * (y46-y15) - 1.0;
subject to d46_16:
	0 >= (x46-x16) * (x46-x16) + (y46-y16) * (y46-y16) - 1.0;
subject to d46_17:
	0 >= (x46-x17) * (x46-x17) + (y46-y17) * (y46-y17) - 1.0;
subject to d46_18:
	0 >= (x46-x18) * (x46-x18) + (y46-y18) * (y46-y18) - 1.0;
subject to d46_19:
	0 >= (x46-x19) * (x46-x19) + (y46-y19) * (y46-y19) - 1.0;
subject to d46_20:
	0 >= (x46-x20) * (x46-x20) + (y46-y20) * (y46-y20) - 1.0;
subject to d46_21:
	0 >= (x46-x21) * (x46-x21) + (y46-y21) * (y46-y21) - 1.0;
subject to d46_22:
	0 >= (x46-x22) * (x46-x22) + (y46-y22) * (y46-y22) - 1.0;
subject to d46_23:
	0 >= (x46-x23) * (x46-x23) + (y46-y23) * (y46-y23) - 1.0;
subject to d46_24:
	0 >= (x46-x24) * (x46-x24) + (y46-y24) * (y46-y24) - 1.0;
subject to d46_25:
	0 >= (x46-x25) * (x46-x25) + (y46-y25) * (y46-y25) - 1.0;
subject to d46_26:
	0 >= (x46-x26) * (x46-x26) + (y46-y26) * (y46-y26) - 1.0;
subject to d46_27:
	0 >= (x46-x27) * (x46-x27) + (y46-y27) * (y46-y27) - 1.0;
subject to d46_28:
	0 >= (x46-x28) * (x46-x28) + (y46-y28) * (y46-y28) - 1.0;
subject to d46_29:
	0 >= (x46-x29) * (x46-x29) + (y46-y29) * (y46-y29) - 1.0;
subject to d46_30:
	0 >= (x46-x30) * (x46-x30) + (y46-y30) * (y46-y30) - 1.0;
subject to d46_31:
	0 >= (x46-x31) * (x46-x31) + (y46-y31) * (y46-y31) - 1.0;
subject to d46_32:
	0 >= (x46-x32) * (x46-x32) + (y46-y32) * (y46-y32) - 1.0;
subject to d46_33:
	0 >= (x46-x33) * (x46-x33) + (y46-y33) * (y46-y33) - 1.0;
subject to d46_34:
	0 >= (x46-x34) * (x46-x34) + (y46-y34) * (y46-y34) - 1.0;
subject to d46_35:
	0 >= (x46-x35) * (x46-x35) + (y46-y35) * (y46-y35) - 1.0;
subject to d46_36:
	0 >= (x46-x36) * (x46-x36) + (y46-y36) * (y46-y36) - 1.0;
subject to d46_37:
	0 >= (x46-x37) * (x46-x37) + (y46-y37) * (y46-y37) - 1.0;
subject to d46_38:
	0 >= (x46-x38) * (x46-x38) + (y46-y38) * (y46-y38) - 1.0;
subject to d46_39:
	0 >= (x46-x39) * (x46-x39) + (y46-y39) * (y46-y39) - 1.0;
subject to d46_40:
	0 >= (x46-x40) * (x46-x40) + (y46-y40) * (y46-y40) - 1.0;
subject to d46_41:
	0 >= (x46-x41) * (x46-x41) + (y46-y41) * (y46-y41) - 1.0;
subject to d46_42:
	0 >= (x46-x42) * (x46-x42) + (y46-y42) * (y46-y42) - 1.0;
subject to d46_43:
	0 >= (x46-x43) * (x46-x43) + (y46-y43) * (y46-y43) - 1.0;
subject to d46_44:
	0 >= (x46-x44) * (x46-x44) + (y46-y44) * (y46-y44) - 1.0;
subject to d46_45:
	0 >= (x46-x45) * (x46-x45) + (y46-y45) * (y46-y45) - 1.0;
subject to d47_0:
	0 >= (x47-x0) * (x47-x0) + (y47-y0) * (y47-y0) - 1.0;
subject to d47_1:
	0 >= (x47-x1) * (x47-x1) + (y47-y1) * (y47-y1) - 1.0;
subject to d47_2:
	0 >= (x47-x2) * (x47-x2) + (y47-y2) * (y47-y2) - 1.0;
subject to d47_3:
	0 >= (x47-x3) * (x47-x3) + (y47-y3) * (y47-y3) - 1.0;
subject to d47_4:
	0 >= (x47-x4) * (x47-x4) + (y47-y4) * (y47-y4) - 1.0;
subject to d47_5:
	0 >= (x47-x5) * (x47-x5) + (y47-y5) * (y47-y5) - 1.0;
subject to d47_6:
	0 >= (x47-x6) * (x47-x6) + (y47-y6) * (y47-y6) - 1.0;
subject to d47_7:
	0 >= (x47-x7) * (x47-x7) + (y47-y7) * (y47-y7) - 1.0;
subject to d47_8:
	0 >= (x47-x8) * (x47-x8) + (y47-y8) * (y47-y8) - 1.0;
subject to d47_9:
	0 >= (x47-x9) * (x47-x9) + (y47-y9) * (y47-y9) - 1.0;
subject to d47_10:
	0 >= (x47-x10) * (x47-x10) + (y47-y10) * (y47-y10) - 1.0;
subject to d47_11:
	0 >= (x47-x11) * (x47-x11) + (y47-y11) * (y47-y11) - 1.0;
subject to d47_12:
	0 >= (x47-x12) * (x47-x12) + (y47-y12) * (y47-y12) - 1.0;
subject to d47_13:
	0 >= (x47-x13) * (x47-x13) + (y47-y13) * (y47-y13) - 1.0;
subject to d47_14:
	0 >= (x47-x14) * (x47-x14) + (y47-y14) * (y47-y14) - 1.0;
subject to d47_15:
	0 >= (x47-x15) * (x47-x15) + (y47-y15) * (y47-y15) - 1.0;
subject to d47_16:
	0 >= (x47-x16) * (x47-x16) + (y47-y16) * (y47-y16) - 1.0;
subject to d47_17:
	0 >= (x47-x17) * (x47-x17) + (y47-y17) * (y47-y17) - 1.0;
subject to d47_18:
	0 >= (x47-x18) * (x47-x18) + (y47-y18) * (y47-y18) - 1.0;
subject to d47_19:
	0 >= (x47-x19) * (x47-x19) + (y47-y19) * (y47-y19) - 1.0;
subject to d47_20:
	0 >= (x47-x20) * (x47-x20) + (y47-y20) * (y47-y20) - 1.0;
subject to d47_21:
	0 >= (x47-x21) * (x47-x21) + (y47-y21) * (y47-y21) - 1.0;
subject to d47_22:
	0 >= (x47-x22) * (x47-x22) + (y47-y22) * (y47-y22) - 1.0;
subject to d47_23:
	0 >= (x47-x23) * (x47-x23) + (y47-y23) * (y47-y23) - 1.0;
subject to d47_24:
	0 >= (x47-x24) * (x47-x24) + (y47-y24) * (y47-y24) - 1.0;
subject to d47_25:
	0 >= (x47-x25) * (x47-x25) + (y47-y25) * (y47-y25) - 1.0;
subject to d47_26:
	0 >= (x47-x26) * (x47-x26) + (y47-y26) * (y47-y26) - 1.0;
subject to d47_27:
	0 >= (x47-x27) * (x47-x27) + (y47-y27) * (y47-y27) - 1.0;
subject to d47_28:
	0 >= (x47-x28) * (x47-x28) + (y47-y28) * (y47-y28) - 1.0;
subject to d47_29:
	0 >= (x47-x29) * (x47-x29) + (y47-y29) * (y47-y29) - 1.0;
subject to d47_30:
	0 >= (x47-x30) * (x47-x30) + (y47-y30) * (y47-y30) - 1.0;
subject to d47_31:
	0 >= (x47-x31) * (x47-x31) + (y47-y31) * (y47-y31) - 1.0;
subject to d47_32:
	0 >= (x47-x32) * (x47-x32) + (y47-y32) * (y47-y32) - 1.0;
subject to d47_33:
	0 >= (x47-x33) * (x47-x33) + (y47-y33) * (y47-y33) - 1.0;
subject to d47_34:
	0 >= (x47-x34) * (x47-x34) + (y47-y34) * (y47-y34) - 1.0;
subject to d47_35:
	0 >= (x47-x35) * (x47-x35) + (y47-y35) * (y47-y35) - 1.0;
subject to d47_36:
	0 >= (x47-x36) * (x47-x36) + (y47-y36) * (y47-y36) - 1.0;
subject to d47_37:
	0 >= (x47-x37) * (x47-x37) + (y47-y37) * (y47-y37) - 1.0;
subject to d47_38:
	0 >= (x47-x38) * (x47-x38) + (y47-y38) * (y47-y38) - 1.0;
subject to d47_39:
	0 >= (x47-x39) * (x47-x39) + (y47-y39) * (y47-y39) - 1.0;
subject to d47_40:
	0 >= (x47-x40) * (x47-x40) + (y47-y40) * (y47-y40) - 1.0;
subject to d47_41:
	0 >= (x47-x41) * (x47-x41) + (y47-y41) * (y47-y41) - 1.0;
subject to d47_42:
	0 >= (x47-x42) * (x47-x42) + (y47-y42) * (y47-y42) - 1.0;
subject to d47_43:
	0 >= (x47-x43) * (x47-x43) + (y47-y43) * (y47-y43) - 1.0;
subject to d47_44:
	0 >= (x47-x44) * (x47-x44) + (y47-y44) * (y47-y44) - 1.0;
subject to d47_45:
	0 >= (x47-x45) * (x47-x45) + (y47-y45) * (y47-y45) - 1.0;
subject to d47_46:
	0 >= (x47-x46) * (x47-x46) + (y47-y46) * (y47-y46) - 1.0;
subject to d48_0:
	0 >= (x48-x0) * (x48-x0) + (y48-y0) * (y48-y0) - 1.0;
subject to d48_1:
	0 >= (x48-x1) * (x48-x1) + (y48-y1) * (y48-y1) - 1.0;
subject to d48_2:
	0 >= (x48-x2) * (x48-x2) + (y48-y2) * (y48-y2) - 1.0;
subject to d48_3:
	0 >= (x48-x3) * (x48-x3) + (y48-y3) * (y48-y3) - 1.0;
subject to d48_4:
	0 >= (x48-x4) * (x48-x4) + (y48-y4) * (y48-y4) - 1.0;
subject to d48_5:
	0 >= (x48-x5) * (x48-x5) + (y48-y5) * (y48-y5) - 1.0;
subject to d48_6:
	0 >= (x48-x6) * (x48-x6) + (y48-y6) * (y48-y6) - 1.0;
subject to d48_7:
	0 >= (x48-x7) * (x48-x7) + (y48-y7) * (y48-y7) - 1.0;
subject to d48_8:
	0 >= (x48-x8) * (x48-x8) + (y48-y8) * (y48-y8) - 1.0;
subject to d48_9:
	0 >= (x48-x9) * (x48-x9) + (y48-y9) * (y48-y9) - 1.0;
subject to d48_10:
	0 >= (x48-x10) * (x48-x10) + (y48-y10) * (y48-y10) - 1.0;
subject to d48_11:
	0 >= (x48-x11) * (x48-x11) + (y48-y11) * (y48-y11) - 1.0;
subject to d48_12:
	0 >= (x48-x12) * (x48-x12) + (y48-y12) * (y48-y12) - 1.0;
subject to d48_13:
	0 >= (x48-x13) * (x48-x13) + (y48-y13) * (y48-y13) - 1.0;
subject to d48_14:
	0 >= (x48-x14) * (x48-x14) + (y48-y14) * (y48-y14) - 1.0;
subject to d48_15:
	0 >= (x48-x15) * (x48-x15) + (y48-y15) * (y48-y15) - 1.0;
subject to d48_16:
	0 >= (x48-x16) * (x48-x16) + (y48-y16) * (y48-y16) - 1.0;
subject to d48_17:
	0 >= (x48-x17) * (x48-x17) + (y48-y17) * (y48-y17) - 1.0;
subject to d48_18:
	0 >= (x48-x18) * (x48-x18) + (y48-y18) * (y48-y18) - 1.0;
subject to d48_19:
	0 >= (x48-x19) * (x48-x19) + (y48-y19) * (y48-y19) - 1.0;
subject to d48_20:
	0 >= (x48-x20) * (x48-x20) + (y48-y20) * (y48-y20) - 1.0;
subject to d48_21:
	0 >= (x48-x21) * (x48-x21) + (y48-y21) * (y48-y21) - 1.0;
subject to d48_22:
	0 >= (x48-x22) * (x48-x22) + (y48-y22) * (y48-y22) - 1.0;
subject to d48_23:
	0 >= (x48-x23) * (x48-x23) + (y48-y23) * (y48-y23) - 1.0;
subject to d48_24:
	0 >= (x48-x24) * (x48-x24) + (y48-y24) * (y48-y24) - 1.0;
subject to d48_25:
	0 >= (x48-x25) * (x48-x25) + (y48-y25) * (y48-y25) - 1.0;
subject to d48_26:
	0 >= (x48-x26) * (x48-x26) + (y48-y26) * (y48-y26) - 1.0;
subject to d48_27:
	0 >= (x48-x27) * (x48-x27) + (y48-y27) * (y48-y27) - 1.0;
subject to d48_28:
	0 >= (x48-x28) * (x48-x28) + (y48-y28) * (y48-y28) - 1.0;
subject to d48_29:
	0 >= (x48-x29) * (x48-x29) + (y48-y29) * (y48-y29) - 1.0;
subject to d48_30:
	0 >= (x48-x30) * (x48-x30) + (y48-y30) * (y48-y30) - 1.0;
subject to d48_31:
	0 >= (x48-x31) * (x48-x31) + (y48-y31) * (y48-y31) - 1.0;
subject to d48_32:
	0 >= (x48-x32) * (x48-x32) + (y48-y32) * (y48-y32) - 1.0;
subject to d48_33:
	0 >= (x48-x33) * (x48-x33) + (y48-y33) * (y48-y33) - 1.0;
subject to d48_34:
	0 >= (x48-x34) * (x48-x34) + (y48-y34) * (y48-y34) - 1.0;
subject to d48_35:
	0 >= (x48-x35) * (x48-x35) + (y48-y35) * (y48-y35) - 1.0;
subject to d48_36:
	0 >= (x48-x36) * (x48-x36) + (y48-y36) * (y48-y36) - 1.0;
subject to d48_37:
	0 >= (x48-x37) * (x48-x37) + (y48-y37) * (y48-y37) - 1.0;
subject to d48_38:
	0 >= (x48-x38) * (x48-x38) + (y48-y38) * (y48-y38) - 1.0;
subject to d48_39:
	0 >= (x48-x39) * (x48-x39) + (y48-y39) * (y48-y39) - 1.0;
subject to d48_40:
	0 >= (x48-x40) * (x48-x40) + (y48-y40) * (y48-y40) - 1.0;
subject to d48_41:
	0 >= (x48-x41) * (x48-x41) + (y48-y41) * (y48-y41) - 1.0;
subject to d48_42:
	0 >= (x48-x42) * (x48-x42) + (y48-y42) * (y48-y42) - 1.0;
subject to d48_43:
	0 >= (x48-x43) * (x48-x43) + (y48-y43) * (y48-y43) - 1.0;
subject to d48_44:
	0 >= (x48-x44) * (x48-x44) + (y48-y44) * (y48-y44) - 1.0;
subject to d48_45:
	0 >= (x48-x45) * (x48-x45) + (y48-y45) * (y48-y45) - 1.0;
subject to d48_46:
	0 >= (x48-x46) * (x48-x46) + (y48-y46) * (y48-y46) - 1.0;
subject to d48_47:
	0 >= (x48-x47) * (x48-x47) + (y48-y47) * (y48-y47) - 1.0;
subject to d49_0:
	0 >= (x49-x0) * (x49-x0) + (y49-y0) * (y49-y0) - 1.0;
subject to d49_1:
	0 >= (x49-x1) * (x49-x1) + (y49-y1) * (y49-y1) - 1.0;
subject to d49_2:
	0 >= (x49-x2) * (x49-x2) + (y49-y2) * (y49-y2) - 1.0;
subject to d49_3:
	0 >= (x49-x3) * (x49-x3) + (y49-y3) * (y49-y3) - 1.0;
subject to d49_4:
	0 >= (x49-x4) * (x49-x4) + (y49-y4) * (y49-y4) - 1.0;
subject to d49_5:
	0 >= (x49-x5) * (x49-x5) + (y49-y5) * (y49-y5) - 1.0;
subject to d49_6:
	0 >= (x49-x6) * (x49-x6) + (y49-y6) * (y49-y6) - 1.0;
subject to d49_7:
	0 >= (x49-x7) * (x49-x7) + (y49-y7) * (y49-y7) - 1.0;
subject to d49_8:
	0 >= (x49-x8) * (x49-x8) + (y49-y8) * (y49-y8) - 1.0;
subject to d49_9:
	0 >= (x49-x9) * (x49-x9) + (y49-y9) * (y49-y9) - 1.0;
subject to d49_10:
	0 >= (x49-x10) * (x49-x10) + (y49-y10) * (y49-y10) - 1.0;
subject to d49_11:
	0 >= (x49-x11) * (x49-x11) + (y49-y11) * (y49-y11) - 1.0;
subject to d49_12:
	0 >= (x49-x12) * (x49-x12) + (y49-y12) * (y49-y12) - 1.0;
subject to d49_13:
	0 >= (x49-x13) * (x49-x13) + (y49-y13) * (y49-y13) - 1.0;
subject to d49_14:
	0 >= (x49-x14) * (x49-x14) + (y49-y14) * (y49-y14) - 1.0;
subject to d49_15:
	0 >= (x49-x15) * (x49-x15) + (y49-y15) * (y49-y15) - 1.0;
subject to d49_16:
	0 >= (x49-x16) * (x49-x16) + (y49-y16) * (y49-y16) - 1.0;
subject to d49_17:
	0 >= (x49-x17) * (x49-x17) + (y49-y17) * (y49-y17) - 1.0;
subject to d49_18:
	0 >= (x49-x18) * (x49-x18) + (y49-y18) * (y49-y18) - 1.0;
subject to d49_19:
	0 >= (x49-x19) * (x49-x19) + (y49-y19) * (y49-y19) - 1.0;
subject to d49_20:
	0 >= (x49-x20) * (x49-x20) + (y49-y20) * (y49-y20) - 1.0;
subject to d49_21:
	0 >= (x49-x21) * (x49-x21) + (y49-y21) * (y49-y21) - 1.0;
subject to d49_22:
	0 >= (x49-x22) * (x49-x22) + (y49-y22) * (y49-y22) - 1.0;
subject to d49_23:
	0 >= (x49-x23) * (x49-x23) + (y49-y23) * (y49-y23) - 1.0;
subject to d49_24:
	0 >= (x49-x24) * (x49-x24) + (y49-y24) * (y49-y24) - 1.0;
subject to d49_25:
	0 >= (x49-x25) * (x49-x25) + (y49-y25) * (y49-y25) - 1.0;
subject to d49_26:
	0 >= (x49-x26) * (x49-x26) + (y49-y26) * (y49-y26) - 1.0;
subject to d49_27:
	0 >= (x49-x27) * (x49-x27) + (y49-y27) * (y49-y27) - 1.0;
subject to d49_28:
	0 >= (x49-x28) * (x49-x28) + (y49-y28) * (y49-y28) - 1.0;
subject to d49_29:
	0 >= (x49-x29) * (x49-x29) + (y49-y29) * (y49-y29) - 1.0;
subject to d49_30:
	0 >= (x49-x30) * (x49-x30) + (y49-y30) * (y49-y30) - 1.0;
subject to d49_31:
	0 >= (x49-x31) * (x49-x31) + (y49-y31) * (y49-y31) - 1.0;
subject to d49_32:
	0 >= (x49-x32) * (x49-x32) + (y49-y32) * (y49-y32) - 1.0;
subject to d49_33:
	0 >= (x49-x33) * (x49-x33) + (y49-y33) * (y49-y33) - 1.0;
subject to d49_34:
	0 >= (x49-x34) * (x49-x34) + (y49-y34) * (y49-y34) - 1.0;
subject to d49_35:
	0 >= (x49-x35) * (x49-x35) + (y49-y35) * (y49-y35) - 1.0;
subject to d49_36:
	0 >= (x49-x36) * (x49-x36) + (y49-y36) * (y49-y36) - 1.0;
subject to d49_37:
	0 >= (x49-x37) * (x49-x37) + (y49-y37) * (y49-y37) - 1.0;
subject to d49_38:
	0 >= (x49-x38) * (x49-x38) + (y49-y38) * (y49-y38) - 1.0;
subject to d49_39:
	0 >= (x49-x39) * (x49-x39) + (y49-y39) * (y49-y39) - 1.0;
subject to d49_40:
	0 >= (x49-x40) * (x49-x40) + (y49-y40) * (y49-y40) - 1.0;
subject to d49_41:
	0 >= (x49-x41) * (x49-x41) + (y49-y41) * (y49-y41) - 1.0;
subject to d49_42:
	0 >= (x49-x42) * (x49-x42) + (y49-y42) * (y49-y42) - 1.0;
subject to d49_43:
	0 >= (x49-x43) * (x49-x43) + (y49-y43) * (y49-y43) - 1.0;
subject to d49_44:
	0 >= (x49-x44) * (x49-x44) + (y49-y44) * (y49-y44) - 1.0;
subject to d49_45:
	0 >= (x49-x45) * (x49-x45) + (y49-y45) * (y49-y45) - 1.0;
subject to d49_46:
	0 >= (x49-x46) * (x49-x46) + (y49-y46) * (y49-y46) - 1.0;
subject to d49_47:
	0 >= (x49-x47) * (x49-x47) + (y49-y47) * (y49-y47) - 1.0;
subject to d49_48:
	0 >= (x49-x48) * (x49-x48) + (y49-y48) * (y49-y48) - 1.0;

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
	display x11;
	display y11;
	display x12;
	display y12;
	display x13;
	display y13;
	display x14;
	display y14;
	display x15;
	display y15;
	display x16;
	display y16;
	display x17;
	display y17;
	display x18;
	display y18;
	display x19;
	display y19;
	display x20;
	display y20;
	display x21;
	display y21;
	display x22;
	display y22;
	display x23;
	display y23;
	display x24;
	display y24;
	display x25;
	display y25;
	display x26;
	display y26;
	display x27;
	display y27;
	display x28;
	display y28;
	display x29;
	display y29;
	display x30;
	display y30;
	display x31;
	display y31;
	display x32;
	display y32;
	display x33;
	display y33;
	display x34;
	display y34;
	display x35;
	display y35;
	display x36;
	display y36;
	display x37;
	display y37;
	display x38;
	display y38;
	display x39;
	display y39;
	display x40;
	display y40;
	display x41;
	display y41;
	display x42;
	display y42;
	display x43;
	display y43;
	display x44;
	display y44;
	display x45;
	display y45;
	display x46;
	display y46;
	display x47;
	display y47;
	display x48;
	display y48;
	display x49;
	display y49;
display obj;
