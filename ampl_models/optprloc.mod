#
#**************************
# SET UP THE INITIAL DATA *
#**************************
#   Problem :
#   *********
#   Optimal positioning of a new product in a multiattribute space.
#   Consider a market of M existing products, a set of N consumers
#   in a multiattribute (dim K) space.
#   Source: Test problem 4 in M. Duran & I.E. Grossmann,
#   "An outer approximation algorithm for a class of mixed integer nonlinear
#    programs", Mathematical Programming 36, pp. 307-339, 1986.
#   SIF input: S. Leyffer, October 1997
#   classification QQR2-AN-30-30
# ... ideal points Z(I,K)
# ... attribute weights W(I,K)
# ... existing products DEL(J,K)
# ... R(I) = min ( sum_K W(I,K) * ( DEL(J,K) - Z(I,K) )**2 (see AUXIL.SIF)
# ... ellipsoid defining closeness to ideal consumer point
# ... definition of attribute space (linear constraints)
# ... additional bounds (one problem in Branch-and-Bound tree)
#XL OPTPRLOC  Y22       1.0
#XU OPTPRLOC  Y4        0.0
#XU OPTPRLOC  Y23       0.0
	param k := 5;
	param m := 10;
	param n := 25;
	param h := 1000.0;
	param z1_1 := 2.26;
	param z1_2 := 5.15;
	param z1_3 := 4.03;
	param z1_4 := 1.74;
	param z1_5 := 4.74;
	param z2_1 := 5.51;
	param z2_2 := 9.01;
	param z2_3 := 3.84;
	param z2_4 := 1.47;
	param z2_5 := 9.92;
	param z3_1 := 4.06;
	param z3_2 := 1.8;
	param z3_3 := 0.71;
	param z3_4 := 9.09;
	param z3_5 := 8.13;
	param z4_1 := 6.3;
	param z4_2 := 0.11;
	param z4_3 := 4.08;
	param z4_4 := 7.29;
	param z4_5 := 4.24;
	param z5_1 := 2.81;
	param z5_2 := 1.65;
	param z5_3 := 8.08;
	param z5_4 := 3.99;
	param z5_5 := 3.51;
	param z6_1 := 4.29;
	param z6_2 := 9.49;
	param z6_3 := 2.24;
	param z6_4 := 9.78;
	param z6_5 := 1.52;
	param z7_1 := 9.76;
	param z7_2 := 3.64;
	param z7_3 := 6.62;
	param z7_4 := 3.66;
	param z7_5 := 9.08;
	param z8_1 := 1.37;
	param z8_2 := 6.99;
	param z8_3 := 7.19;
	param z8_4 := 3.03;
	param z8_5 := 3.39;
	param z9_1 := 8.89;
	param z9_2 := 8.29;
	param z9_3 := 6.05;
	param z9_4 := 7.48;
	param z9_5 := 4.09;
	param z10_1 := 7.42;
	param z10_2 := 4.6;
	param z10_3 := 0.3;
	param z10_4 := 0.97;
	param z10_5 := 8.77;
	param z11_1 := 1.54;
	param z11_2 := 7.06;
	param z11_3 := 0.01;
	param z11_4 := 1.23;
	param z11_5 := 3.11;
	param z12_1 := 7.74;
	param z12_2 := 4.4;
	param z12_3 := 7.93;
	param z12_4 := 5.95;
	param z12_5 := 4.88;
	param z13_1 := 9.94;
	param z13_2 := 5.21;
	param z13_3 := 8.58;
	param z13_4 := 0.13;
	param z13_5 := 4.57;
	param z14_1 := 9.54;
	param z14_2 := 1.57;
	param z14_3 := 9.66;
	param z14_4 := 5.24;
	param z14_5 := 7.9;
	param z15_1 := 7.46;
	param z15_2 := 8.81;
	param z15_3 := 1.67;
	param z15_4 := 6.47;
	param z15_5 := 1.81;
	param z16_1 := 0.56;
	param z16_2 := 8.1;
	param z16_3 := 0.19;
	param z16_4 := 6.11;
	param z16_5 := 6.4;
	param z17_1 := 3.86;
	param z17_2 := 6.68;
	param z17_3 := 6.42;
	param z17_4 := 7.29;
	param z17_5 := 4.66;
	param z18_1 := 2.98;
	param z18_2 := 2.98;
	param z18_3 := 3.03;
	param z18_4 := 0.02;
	param z18_5 := 0.67;
	param z19_1 := 3.61;
	param z19_2 := 7.62;
	param z19_3 := 1.79;
	param z19_4 := 7.8;
	param z19_5 := 9.81;
	param z20_1 := 5.68;
	param z20_2 := 4.24;
	param z20_3 := 4.17;
	param z20_4 := 6.75;
	param z20_5 := 1.08;
	param z21_1 := 5.48;
	param z21_2 := 3.74;
	param z21_3 := 3.34;
	param z21_4 := 6.22;
	param z21_5 := 7.94;
	param z22_1 := 8.13;
	param z22_2 := 8.72;
	param z22_3 := 3.93;
	param z22_4 := 8.8;
	param z22_5 := 8.56;
	param z23_1 := 1.37;
	param z23_2 := 0.54;
	param z23_3 := 1.55;
	param z23_4 := 5.56;
	param z23_5 := 5.85;
	param z24_1 := 8.79;
	param z24_2 := 5.04;
	param z24_3 := 4.83;
	param z24_4 := 6.94;
	param z24_5 := 0.38;
	param z25_1 := 2.66;
	param z25_2 := 4.19;
	param z25_3 := 6.49;
	param z25_4 := 8.04;
	param z25_5 := 1.66;
	param w1_1 := 9.57;
	param w1_2 := 2.74;
	param w1_3 := 9.75;
	param w1_4 := 3.96;
	param w1_5 := 8.67;
	param w2_1 := 8.38;
	param w2_2 := 3.93;
	param w2_3 := 5.18;
	param w2_4 := 5.2;
	param w2_5 := 7.82;
	param w3_1 := 9.81;
	param w3_2 := 0.04;
	param w3_3 := 4.21;
	param w3_4 := 7.38;
	param w3_5 := 4.11;
	param w4_1 := 7.41;
	param w4_2 := 6.08;
	param w4_3 := 5.46;
	param w4_4 := 4.86;
	param w4_5 := 1.48;
	param w5_1 := 9.96;
	param w5_2 := 9.13;
	param w5_3 := 2.95;
	param w5_4 := 8.25;
	param w5_5 := 3.58;
	param w6_1 := 9.39;
	param w6_2 := 4.27;
	param w6_3 := 5.09;
	param w6_4 := 1.81;
	param w6_5 := 7.58;
	param w7_1 := 1.88;
	param w7_2 := 7.2;
	param w7_3 := 6.65;
	param w7_4 := 1.74;
	param w7_5 := 2.86;
	param w8_1 := 4.01;
	param w8_2 := 2.67;
	param w8_3 := 4.86;
	param w8_4 := 2.55;
	param w8_5 := 6.91;
	param w9_1 := 4.18;
	param w9_2 := 1.92;
	param w9_3 := 2.6;
	param w9_4 := 7.15;
	param w9_5 := 2.86;
	param w10_1 := 7.81;
	param w10_2 := 2.14;
	param w10_3 := 9.63;
	param w10_4 := 7.61;
	param w10_5 := 9.17;
	param w11_1 := 8.96;
	param w11_2 := 3.47;
	param w11_3 := 5.49;
	param w11_4 := 4.73;
	param w11_5 := 9.43;
	param w12_1 := 9.94;
	param w12_2 := 1.63;
	param w12_3 := 1.23;
	param w12_4 := 4.33;
	param w12_5 := 7.08;
	param w13_1 := 0.31;
	param w13_2 := 5.0;
	param w13_3 := 0.16;
	param w13_4 := 2.52;
	param w13_5 := 3.08;
	param w14_1 := 6.02;
	param w14_2 := 0.92;
	param w14_3 := 7.47;
	param w14_4 := 9.74;
	param w14_5 := 1.76;
	param w15_1 := 5.06;
	param w15_2 := 4.52;
	param w15_3 := 1.89;
	param w15_4 := 1.22;
	param w15_5 := 9.05;
	param w16_1 := 5.92;
	param w16_2 := 2.56;
	param w16_3 := 7.74;
	param w16_4 := 6.96;
	param w16_5 := 5.18;
	param w17_1 := 6.45;
	param w17_2 := 1.52;
	param w17_3 := 0.06;
	param w17_4 := 5.34;
	param w17_5 := 8.47;
	param w18_1 := 1.04;
	param w18_2 := 1.36;
	param w18_3 := 5.99;
	param w18_4 := 8.1;
	param w18_5 := 5.22;
	param w19_1 := 1.4;
	param w19_2 := 1.35;
	param w19_3 := 0.59;
	param w19_4 := 8.58;
	param w19_5 := 1.21;
	param w20_1 := 6.68;
	param w20_2 := 9.48;
	param w20_3 := 1.6;
	param w20_4 := 6.74;
	param w20_5 := 8.92;
	param w21_1 := 1.95;
	param w21_2 := 0.46;
	param w21_3 := 2.9;
	param w21_4 := 1.79;
	param w21_5 := 0.99;
	param w22_1 := 5.18;
	param w22_2 := 5.1;
	param w22_3 := 8.81;
	param w22_4 := 3.27;
	param w22_5 := 9.63;
	param w23_1 := 1.47;
	param w23_2 := 5.71;
	param w23_3 := 6.95;
	param w23_4 := 1.42;
	param w23_5 := 3.49;
	param w24_1 := 5.4;
	param w24_2 := 3.12;
	param w24_3 := 5.37;
	param w24_4 := 6.1;
	param w24_5 := 3.71;
	param w25_1 := 6.32;
	param w25_2 := 0.81;
	param w25_3 := 6.12;
	param w25_4 := 6.73;
	param w25_5 := 7.93;
	param del1_1 := 0.62;
	param del1_2 := 5.06;
	param del1_3 := 7.82;
	param del1_4 := 0.22;
	param del1_5 := 4.42;
	param del2_1 := 5.21;
	param del2_2 := 2.66;
	param del2_3 := 9.54;
	param del2_4 := 5.03;
	param del2_5 := 8.01;
	param del3_1 := 5.27;
	param del3_2 := 7.72;
	param del3_3 := 7.97;
	param del3_4 := 3.31;
	param del3_5 := 6.56;
	param del4_1 := 1.02;
	param del4_2 := 8.89;
	param del4_3 := 8.77;
	param del4_4 := 3.1;
	param del4_5 := 6.66;
	param del5_1 := 1.26;
	param del5_2 := 6.8;
	param del5_3 := 2.3;
	param del5_4 := 1.75;
	param del5_5 := 6.65;
	param del6_1 := 3.74;
	param del6_2 := 9.06;
	param del6_3 := 9.8;
	param del6_4 := 3.01;
	param del6_5 := 9.52;
	param del7_1 := 4.64;
	param del7_2 := 7.99;
	param del7_3 := 6.69;
	param del7_4 := 5.88;
	param del7_5 := 8.23;
	param del8_1 := 8.35;
	param del8_2 := 3.79;
	param del8_3 := 1.19;
	param del8_4 := 1.96;
	param del8_5 := 5.88;
	param del9_1 := 6.44;
	param del9_2 := 0.17;
	param del9_3 := 9.93;
	param del9_4 := 6.8;
	param del9_5 := 9.75;
	param del10_1 := 6.49;
	param del10_2 := 1.92;
	param del10_3 := 0.05;
	param del10_4 := 4.89;
	param del10_5 := 6.43;
	param r1 := 77.83985;
	param r2 := 175.971;
	param r3 := 201.8226;
	param r4 := 143.9533;
	param r5 := 154.3895;
	param r6 := 433.3177;
	param r7 := 109.0764;
	param r8 := 41.59592;
	param r9 := 144.0623;
	param r10 := 99.83416;
	param r11 := 149.1791;
	param r12 := 123.8074;
	param r13 := 27.22197;
	param r14 := 89.92683;
	param r15 := 293.0766;
	param r16 := 174.317;
	param r17 := 125.1028;
	param r18 := 222.8417;
	param r19 := 50.48593;
	param r20 := 361.1973;
	param r21 := 40.32642;
	param r22 := 161.8518;
	param r23 := 66.85827;
	param r24 := 340.5807;
	param r25 := 407.52;
	param rph1 := (77.83985) + (1000.0);
	param rph2 := (175.971) + (1000.0);
	param rph3 := (201.8226) + (1000.0);
	param rph4 := (143.9533) + (1000.0);
	param rph5 := (154.3895) + (1000.0);
	param rph6 := (433.3177) + (1000.0);
	param rph7 := (109.0764) + (1000.0);
	param rph8 := (41.59592) + (1000.0);
	param rph9 := (144.0623) + (1000.0);
	param rph10 := (99.83416) + (1000.0);
	param rph11 := (149.1791) + (1000.0);
	param rph12 := (123.8074) + (1000.0);
	param rph13 := (27.22197) + (1000.0);
	param rph14 := (89.92683) + (1000.0);
	param rph15 := (293.0766) + (1000.0);
	param rph16 := (174.317) + (1000.0);
	param rph17 := (125.1028) + (1000.0);
	param rph18 := (222.8417) + (1000.0);
	param rph19 := (50.48593) + (1000.0);
	param rph20 := (361.1973) + (1000.0);
	param rph21 := (40.32642) + (1000.0);
	param rph22 := (161.8518) + (1000.0);
	param rph23 := (66.85827) + (1000.0);
	param rph24 := (340.5807) + (1000.0);
	param rph25 := (407.52) + (1000.0);

	var x1 >= 2.0 ,  <= 4.5, := 0;
	var x2 >= 0.0 ,  <= 8.0, := 0;
	var x3 >= 3.0 ,  <= 9.0, := 0;
	var x4 >= 0.0 ,  <= 5.0, := 0;
	var x5 >= 4.0 ,  <= 10.0, := 0;
	var y1 >= 0.0 ,  <= 1.0, := 0;
	var y2 >= 0.0 ,  <= 1.0, := 0;
	var y3 >= 0.0 ,  <= 1.0, := 0;
	var y4 >= 0.0 ,  <= 1.0, := 0;
	var y5 >= 0.0 ,  <= 1.0, := 0;
	var y6 >= 0.0 ,  <= 1.0, := 0;
	var y7 >= 0.0 ,  <= 1.0, := 0;
	var y8 >= 0.0 ,  <= 1.0, := 0;
	var y9 >= 0.0 ,  <= 1.0, := 0;
	var y10 >= 0.0 ,  <= 1.0, := 0;
	var y11 >= 0.0 ,  <= 1.0, := 0;
	var y12 >= 0.0 ,  <= 1.0, := 0;
	var y13 >= 0.0 ,  <= 1.0, := 0;
	var y14 >= 0.0 ,  <= 1.0, := 0;
	var y15 >= 0.0 ,  <= 1.0, := 0;
	var y16 >= 0.0 ,  <= 1.0, := 0;
	var y17 >= 0.0 ,  <= 1.0, := 0;
	var y18 >= 0.0 ,  <= 1.0, := 0;
	var y19 >= 0.0 ,  <= 1.0, := 0;
	var y20 >= 0.0 ,  <= 1.0, := 0;
	var y21 >= 0.0 ,  <= 1.0, := 0;
	var y22 >= 0.0 ,  <= 1.0, := 0;
	var y23 >= 0.0 ,  <= 1.0, := 0;
	var y24 >= 0.0 ,  <= 1.0, := 0;
	var y25 >= 0.0 ,  <= 1.0, := 0;

minimize obj:
	0.6*(x1 - 0.0 ) ^2 + 0.1*(x4 - 0.0 ) ^2 - y1 - 0.2*y2 - y3 - 0.2*y4 - 0.9*y5 - 
	0.9*y6 - 0.1*y7 - 0.8*y8 - y9 - 0.4*y10 - y11 - 0.3*y12 - 0.1*y13 - 0.3*y14 - 
	0.5*y15 - 0.9*y16 - 0.8*y17 - 0.1*y18 - 0.9*y19 - y20 - y21 - y22 - 0.2*y23 - 
	0.7*y24 - 0.7*y25 - 0.9*x2 - 0.5*x3 + x5;

subject to elli1:
	0 >= 9.57*(x1 - 2.26 ) ^2 + 2.74*(x2 - 5.15 ) ^2 + 9.75*(x3 - 4.03 ) ^2 + 
	3.96*(x4 - 1.74 ) ^2 + 8.67*(x5 - 4.74 ) ^2 + 1000.0*y1 - 1077.83985;
subject to elli2:
	0 >= 8.38*(x1 - 5.51 ) ^2 + 3.93*(x2 - 9.01 ) ^2 + 5.18*(x3 - 3.84 ) ^2 + 
	5.2*(x4 - 1.47 ) ^2 + 7.82*(x5 - 9.92 ) ^2 + 1000.0*y2 - 1175.971;
subject to elli3:
	0 >= 9.81*(x1 - 4.06 ) ^2 + 0.04*(x2 - 1.8 ) ^2 + 4.21*(x3 - 0.71 ) ^2 + 
	7.38*(x4 - 9.09 ) ^2 + 4.11*(x5 - 8.13 ) ^2 + 1000.0*y3 - 1201.8226;
subject to elli4:
	0 >= 7.41*(x1 - 6.3 ) ^2 + 6.08*(x2 - 0.11 ) ^2 + 5.46*(x3 - 4.08 ) ^2 + 
	4.86*(x4 - 7.29 ) ^2 + 1.48*(x5 - 4.24 ) ^2 + 1000.0*y4 - 1143.9533000000001;
subject to elli5:
	0 >= 9.96*(x1 - 2.81 ) ^2 + 9.13*(x2 - 1.65 ) ^2 + 2.95*(x3 - 8.08 ) ^2 + 
	8.25*(x4 - 3.99 ) ^2 + 3.58*(x5 - 3.51 ) ^2 + 1000.0*y5 - 1154.3895;
subject to elli6:
	0 >= 9.39*(x1 - 4.29 ) ^2 + 4.27*(x2 - 9.49 ) ^2 + 5.09*(x3 - 2.24 ) ^2 + 
	1.81*(x4 - 9.78 ) ^2 + 7.58*(x5 - 1.52 ) ^2 + 1000.0*y6 - 1433.3177;
subject to elli7:
	0 >= 1.88*(x1 - 9.76 ) ^2 + 7.2*(x2 - 3.64 ) ^2 + 6.65*(x3 - 6.62 ) ^2 + 
	1.74*(x4 - 3.66 ) ^2 + 2.86*(x5 - 9.08 ) ^2 + 1000.0*y7 - 1109.0764;
subject to elli8:
	0 >= 4.01*(x1 - 1.37 ) ^2 + 2.67*(x2 - 6.99 ) ^2 + 4.86*(x3 - 7.19 ) ^2 + 
	2.55*(x4 - 3.03 ) ^2 + 6.91*(x5 - 3.39 ) ^2 + 1000.0*y8 - 1041.59592;
subject to elli9:
	0 >= 4.18*(x1 - 8.89 ) ^2 + 1.92*(x2 - 8.29 ) ^2 + 2.6*(x3 - 6.05 ) ^2 + 
	7.15*(x4 - 7.48 ) ^2 + 2.86*(x5 - 4.09 ) ^2 + 1000.0*y9 - 1144.0623;
subject to elli10:
	0 >= 7.81*(x1 - 7.42 ) ^2 + 2.14*(x2 - 4.6 ) ^2 + 9.63*(x3 - 0.3 ) ^2 + 
	7.61*(x4 - 0.97 ) ^2 + 9.17*(x5 - 8.77 ) ^2 + 1000.0*y10 - 1099.8341599999999;
subject to elli11:
	0 >= 8.96*(x1 - 1.54 ) ^2 + 3.47*(x2 - 7.06 ) ^2 + 5.49*(x3 - 0.01 ) ^2 + 
	4.73*(x4 - 1.23 ) ^2 + 9.43*(x5 - 3.11 ) ^2 + 1000.0*y11 - 1149.1791;
subject to elli12:
	0 >= 9.94*(x1 - 7.74 ) ^2 + 1.63*(x2 - 4.4 ) ^2 + 1.23*(x3 - 7.93 ) ^2 + 
	4.33*(x4 - 5.95 ) ^2 + 7.08*(x5 - 4.88 ) ^2 + 1000.0*y12 - 1123.8074;
subject to elli13:
	0 >= 0.31*(x1 - 9.94 ) ^2 + 5.0*(x2 - 5.21 ) ^2 + 0.16*(x3 - 8.58 ) ^2 + 
	2.52*(x4 - 0.13 ) ^2 + 3.08*(x5 - 4.57 ) ^2 + 1000.0*y13 - 1027.22197;
subject to elli14:
	0 >= 6.02*(x1 - 9.54 ) ^2 + 0.92*(x2 - 1.57 ) ^2 + 7.47*(x3 - 9.66 ) ^2 + 
	9.74*(x4 - 5.24 ) ^2 + 1.76*(x5 - 7.9 ) ^2 + 1000.0*y14 - 1089.9268299999999;
subject to elli15:
	0 >= 5.06*(x1 - 7.46 ) ^2 + 4.52*(x2 - 8.81 ) ^2 + 1.89*(x3 - 1.67 ) ^2 + 
	1.22*(x4 - 6.47 ) ^2 + 9.05*(x5 - 1.81 ) ^2 + 1000.0*y15 - 1293.0765999999999;
subject to elli16:
	0 >= 5.92*(x1 - 0.56 ) ^2 + 2.56*(x2 - 8.1 ) ^2 + 7.74*(x3 - 0.19 ) ^2 + 
	6.96*(x4 - 6.11 ) ^2 + 5.18*(x5 - 6.4 ) ^2 + 1000.0*y16 - 1174.317;
subject to elli17:
	0 >= 6.45*(x1 - 3.86 ) ^2 + 1.52*(x2 - 6.68 ) ^2 + 0.06*(x3 - 6.42 ) ^2 + 
	5.34*(x4 - 7.29 ) ^2 + 8.47*(x5 - 4.66 ) ^2 + 1000.0*y17 - 1125.1028000000001;
subject to elli18:
	0 >= 1.04*(x1 - 2.98 ) ^2 + 1.36*(x2 - 2.98 ) ^2 + 5.99*(x3 - 3.03 ) ^2 + 
	8.1*(x4 - 0.02 ) ^2 + 5.22*(x5 - 0.67 ) ^2 + 1000.0*y18 - 1222.8417;
subject to elli19:
	0 >= 1.4*(x1 - 3.61 ) ^2 + 1.35*(x2 - 7.62 ) ^2 + 0.59*(x3 - 1.79 ) ^2 + 
	8.58*(x4 - 7.8 ) ^2 + 1.21*(x5 - 9.81 ) ^2 + 1000.0*y19 - 1050.48593;
subject to elli20:
	0 >= 6.68*(x1 - 5.68 ) ^2 + 9.48*(x2 - 4.24 ) ^2 + 1.6*(x3 - 4.17 ) ^2 + 
	6.74*(x4 - 6.75 ) ^2 + 8.92*(x5 - 1.08 ) ^2 + 1000.0*y20 - 1361.1973;
subject to elli21:
	0 >= 1.95*(x1 - 5.48 ) ^2 + 0.46*(x2 - 3.74 ) ^2 + 2.9*(x3 - 3.34 ) ^2 + 
	1.79*(x4 - 6.22 ) ^2 + 0.99*(x5 - 7.94 ) ^2 + 1000.0*y21 - 1040.32642;
subject to elli22:
	0 >= 5.18*(x1 - 8.13 ) ^2 + 5.1*(x2 - 8.72 ) ^2 + 8.81*(x3 - 3.93 ) ^2 + 
	3.27*(x4 - 8.8 ) ^2 + 9.63*(x5 - 8.56 ) ^2 + 1000.0*y22 - 1161.8518;
subject to elli23:
	0 >= 1.47*(x1 - 1.37 ) ^2 + 5.71*(x2 - 0.54 ) ^2 + 6.95*(x3 - 1.55 ) ^2 + 
	1.42*(x4 - 5.56 ) ^2 + 3.49*(x5 - 5.85 ) ^2 + 1000.0*y23 - 1066.85827;
subject to elli24:
	0 >= 5.4*(x1 - 8.79 ) ^2 + 3.12*(x2 - 5.04 ) ^2 + 5.37*(x3 - 4.83 ) ^2 + 
	6.1*(x4 - 6.94 ) ^2 + 3.71*(x5 - 0.38 ) ^2 + 1000.0*y24 - 1340.5807;
subject to elli25:
	0 >= 6.32*(x1 - 2.66 ) ^2 + 0.81*(x2 - 4.19 ) ^2 + 6.12*(x3 - 6.49 ) ^2 + 
	6.73*(x4 - 8.04 ) ^2 + 7.93*(x5 - 1.66 ) ^2 + 1000.0*y25 - 1407.52;
subject to lin1:
	0 >= x1 - x2 + x3 + x4 + x5 - 10.0;
subject to lin2:
	0 >= 0.6*x1 - 0.9*x2 - 0.5*x3 + 0.1*x4 + x5 + 0.64;
subject to lin3:
	0 <= x1 - x2 + x3 - x4 + x5 - 0.69;
subject to lin4:
	0 >= 0.157*x1 + 0.05*x2 - 1.5;
subject to lin5:
	0 <= 0.25*x2 + 1.05*x4 - 0.3*x5 - 4.5;

solve;
	display x1;
	display x2;
	display x3;
	display x4;
	display x5;
	display y1;
	display y2;
	display y3;
	display y4;
	display y5;
	display y6;
	display y7;
	display y8;
	display y9;
	display y10;
	display y11;
	display y12;
	display y13;
	display y14;
	display y15;
	display y16;
	display y17;
	display y18;
	display y19;
	display y20;
	display y21;
	display y22;
	display y23;
	display y24;
	display y25;
display obj;
