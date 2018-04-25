#
#**************************
# SET UP THE INITIAL DATA *
#**************************
#   Problem :
#   *********
#   A  two-dimensional base  flow  problem in an inclined enclosure.
#   Heat flux constant at y = +/- 1
#   Low Reynold's number
#   The flow is considered in a square of length 2,  centered on the
#   origin and aligned with the x-y axes. The square is divided into
#   4 n ** 2  sub-squares,  each of  length 1 / n.  The differential
#   equation is replaced by  discrete nonlinear equations at each of 
#   the grid points. 
#   The differential equation relates the vorticity, temperature and
#   a stream function.
#   
#   Source: 
#   J. N. Shadid
#   "Experimental and computational study of the stability
#   of Natural convection flow in an inclined enclosure",
#   Ph. D. Thesis, University of Minnesota, 1989,
#   problem SP2 (pp.128-130), 
#   SIF input: Nick Gould, August 1993.
#   classification NQR2-MY-V-V
#   Half the number of discretization intervals
#   Number of variables = 3(2M+1)**2 
#IE M                   1              $ n=27
#IE M                   2              $ n=75
#IE M                   5              $ n=363
#IE M                   8              $ n=867
#   Define the Rayleigh number. NB: This determines the difficulty 
#   of the problem.
#   Set pi.
#   Define other problem parameters
# Case 2. Constant heat flux at y = +/- 1
#   Define a few helpful parameters
#   Define a vorticity(OM), temperature(PH) and stream function(PS)
#   variable per discretized point in the square
#   Define three equations per interior node 
#   The stream function equation(S) - linear (6.57a in the thesis).
#   The vorticity equation(V) - linear (6.57b in the thesis).
#   The thermal energy equation(E) - quadratic (6.57c in the thesis).
#  Boundary conditions on the temperature.
#  Boundary conditions on the vorticity. NB: Steady state assumed
	param m := 10;
	param ra := 1000.0;
	param pid4 := atan(1.0);
	param pi := 4.0 * (atan(1.0));
	param ax := 1.0;
	param theta := 0.5 * (4.0 * (atan(1.0)));
	param a1 := 1.0;
	param a2 := 0.0;
	param a3 := -1.0;
	param b1 := 1.0;
	param b2 := 0.0;
	param b3 := -1.0;
	param f1 := 1.0;
	param f2 := 0.0;
	param f3 := 0.0;
	param g1 := 1.0;
	param g2 := 0.0;
	param g3 := 0.0;
	param mm1 := -1 + (10);
	param h := 1.0 / (10.0);
	param h2 := (1.0 / (10.0)) * (1.0 / (10.0));
	param axx := (1.0) * (1.0);
	param sintheta := sin(0.5 * (4.0 * (atan(1.0))));
	param costheta := cos(0.5 * (4.0 * (atan(1.0))));
	param pi1 := -0.5 * (((1.0) * (1000.0)) * (cos(0.5 * (4.0 * (atan(1.0))))));
	param pi2 := 0.5 * ((((1.0) * (1.0)) * (1000.0)) * (sin(0.5 * (4.0 * 
	(atan(1.0))))));
	param axd2 := 0.5 * (1.0);
	param axxd2 := 0.5 * ((1.0) * (1.0));
	param axxd4 := 0.25 * ((1.0) * (1.0));
	param pi1d2h := (-0.5 * (((1.0) * (1000.0)) * (cos(0.5 * (4.0 * 
	(atan(1.0))))))) * (0.5 * (10.0));
	param pi2d2h := (0.5 * ((((1.0) * (1.0)) * (1000.0)) * (sin(0.5 * (4.0 * 
	(atan(1.0))))))) * (0.5 * (10.0));
	param axdh2 := (1.0) * ((10.0) * (10.0));
	param axd4h2 := 0.25 * ((1.0) * ((10.0) * (10.0)));
	param axxdh2 := ((1.0) * (1.0)) * ((10.0) * (10.0));
	param jp := 1 + (9);
	param jm := -1 + (9);
	param ip := 1 + (9);
	param im := -1 + (9);

	var omm10_m10;
	var phm10_m10;
	var psm10_m10 >= 1.0 ,  <= 1.0;
	var omm9_m10;
	var phm9_m10;
	var psm9_m10 >= 1.0 ,  <= 1.0;
	var omm8_m10;
	var phm8_m10;
	var psm8_m10 >= 1.0 ,  <= 1.0;
	var omm7_m10;
	var phm7_m10;
	var psm7_m10 >= 1.0 ,  <= 1.0;
	var omm6_m10;
	var phm6_m10;
	var psm6_m10 >= 1.0 ,  <= 1.0;
	var omm5_m10;
	var phm5_m10;
	var psm5_m10 >= 1.0 ,  <= 1.0;
	var omm4_m10;
	var phm4_m10;
	var psm4_m10 >= 1.0 ,  <= 1.0;
	var omm3_m10;
	var phm3_m10;
	var psm3_m10 >= 1.0 ,  <= 1.0;
	var omm2_m10;
	var phm2_m10;
	var psm2_m10 >= 1.0 ,  <= 1.0;
	var omm1_m10;
	var phm1_m10;
	var psm1_m10 >= 1.0 ,  <= 1.0;
	var om0_m10;
	var ph0_m10;
	var ps0_m10 >= 1.0 ,  <= 1.0;
	var om1_m10;
	var ph1_m10;
	var ps1_m10 >= 1.0 ,  <= 1.0;
	var om2_m10;
	var ph2_m10;
	var ps2_m10 >= 1.0 ,  <= 1.0;
	var om3_m10;
	var ph3_m10;
	var ps3_m10 >= 1.0 ,  <= 1.0;
	var om4_m10;
	var ph4_m10;
	var ps4_m10 >= 1.0 ,  <= 1.0;
	var om5_m10;
	var ph5_m10;
	var ps5_m10 >= 1.0 ,  <= 1.0;
	var om6_m10;
	var ph6_m10;
	var ps6_m10 >= 1.0 ,  <= 1.0;
	var om7_m10;
	var ph7_m10;
	var ps7_m10 >= 1.0 ,  <= 1.0;
	var om8_m10;
	var ph8_m10;
	var ps8_m10 >= 1.0 ,  <= 1.0;
	var om9_m10;
	var ph9_m10;
	var ps9_m10 >= 1.0 ,  <= 1.0;
	var om10_m10;
	var ph10_m10;
	var ps10_m10 >= 1.0 ,  <= 1.0;
	var omm10_m9;
	var phm10_m9;
	var psm10_m9 >= 1.0 ,  <= 1.0;
	var omm9_m9;
	var phm9_m9;
	var psm9_m9;
	var omm8_m9;
	var phm8_m9;
	var psm8_m9;
	var omm7_m9;
	var phm7_m9;
	var psm7_m9;
	var omm6_m9;
	var phm6_m9;
	var psm6_m9;
	var omm5_m9;
	var phm5_m9;
	var psm5_m9;
	var omm4_m9;
	var phm4_m9;
	var psm4_m9;
	var omm3_m9;
	var phm3_m9;
	var psm3_m9;
	var omm2_m9;
	var phm2_m9;
	var psm2_m9;
	var omm1_m9;
	var phm1_m9;
	var psm1_m9;
	var om0_m9;
	var ph0_m9;
	var ps0_m9;
	var om1_m9;
	var ph1_m9;
	var ps1_m9;
	var om2_m9;
	var ph2_m9;
	var ps2_m9;
	var om3_m9;
	var ph3_m9;
	var ps3_m9;
	var om4_m9;
	var ph4_m9;
	var ps4_m9;
	var om5_m9;
	var ph5_m9;
	var ps5_m9;
	var om6_m9;
	var ph6_m9;
	var ps6_m9;
	var om7_m9;
	var ph7_m9;
	var ps7_m9;
	var om8_m9;
	var ph8_m9;
	var ps8_m9;
	var om9_m9;
	var ph9_m9;
	var ps9_m9;
	var om10_m9;
	var ph10_m9;
	var ps10_m9 >= 1.0 ,  <= 1.0;
	var omm10_m8;
	var phm10_m8;
	var psm10_m8 >= 1.0 ,  <= 1.0;
	var omm9_m8;
	var phm9_m8;
	var psm9_m8;
	var omm8_m8;
	var phm8_m8;
	var psm8_m8;
	var omm7_m8;
	var phm7_m8;
	var psm7_m8;
	var omm6_m8;
	var phm6_m8;
	var psm6_m8;
	var omm5_m8;
	var phm5_m8;
	var psm5_m8;
	var omm4_m8;
	var phm4_m8;
	var psm4_m8;
	var omm3_m8;
	var phm3_m8;
	var psm3_m8;
	var omm2_m8;
	var phm2_m8;
	var psm2_m8;
	var omm1_m8;
	var phm1_m8;
	var psm1_m8;
	var om0_m8;
	var ph0_m8;
	var ps0_m8;
	var om1_m8;
	var ph1_m8;
	var ps1_m8;
	var om2_m8;
	var ph2_m8;
	var ps2_m8;
	var om3_m8;
	var ph3_m8;
	var ps3_m8;
	var om4_m8;
	var ph4_m8;
	var ps4_m8;
	var om5_m8;
	var ph5_m8;
	var ps5_m8;
	var om6_m8;
	var ph6_m8;
	var ps6_m8;
	var om7_m8;
	var ph7_m8;
	var ps7_m8;
	var om8_m8;
	var ph8_m8;
	var ps8_m8;
	var om9_m8;
	var ph9_m8;
	var ps9_m8;
	var om10_m8;
	var ph10_m8;
	var ps10_m8 >= 1.0 ,  <= 1.0;
	var omm10_m7;
	var phm10_m7;
	var psm10_m7 >= 1.0 ,  <= 1.0;
	var omm9_m7;
	var phm9_m7;
	var psm9_m7;
	var omm8_m7;
	var phm8_m7;
	var psm8_m7;
	var omm7_m7;
	var phm7_m7;
	var psm7_m7;
	var omm6_m7;
	var phm6_m7;
	var psm6_m7;
	var omm5_m7;
	var phm5_m7;
	var psm5_m7;
	var omm4_m7;
	var phm4_m7;
	var psm4_m7;
	var omm3_m7;
	var phm3_m7;
	var psm3_m7;
	var omm2_m7;
	var phm2_m7;
	var psm2_m7;
	var omm1_m7;
	var phm1_m7;
	var psm1_m7;
	var om0_m7;
	var ph0_m7;
	var ps0_m7;
	var om1_m7;
	var ph1_m7;
	var ps1_m7;
	var om2_m7;
	var ph2_m7;
	var ps2_m7;
	var om3_m7;
	var ph3_m7;
	var ps3_m7;
	var om4_m7;
	var ph4_m7;
	var ps4_m7;
	var om5_m7;
	var ph5_m7;
	var ps5_m7;
	var om6_m7;
	var ph6_m7;
	var ps6_m7;
	var om7_m7;
	var ph7_m7;
	var ps7_m7;
	var om8_m7;
	var ph8_m7;
	var ps8_m7;
	var om9_m7;
	var ph9_m7;
	var ps9_m7;
	var om10_m7;
	var ph10_m7;
	var ps10_m7 >= 1.0 ,  <= 1.0;
	var omm10_m6;
	var phm10_m6;
	var psm10_m6 >= 1.0 ,  <= 1.0;
	var omm9_m6;
	var phm9_m6;
	var psm9_m6;
	var omm8_m6;
	var phm8_m6;
	var psm8_m6;
	var omm7_m6;
	var phm7_m6;
	var psm7_m6;
	var omm6_m6;
	var phm6_m6;
	var psm6_m6;
	var omm5_m6;
	var phm5_m6;
	var psm5_m6;
	var omm4_m6;
	var phm4_m6;
	var psm4_m6;
	var omm3_m6;
	var phm3_m6;
	var psm3_m6;
	var omm2_m6;
	var phm2_m6;
	var psm2_m6;
	var omm1_m6;
	var phm1_m6;
	var psm1_m6;
	var om0_m6;
	var ph0_m6;
	var ps0_m6;
	var om1_m6;
	var ph1_m6;
	var ps1_m6;
	var om2_m6;
	var ph2_m6;
	var ps2_m6;
	var om3_m6;
	var ph3_m6;
	var ps3_m6;
	var om4_m6;
	var ph4_m6;
	var ps4_m6;
	var om5_m6;
	var ph5_m6;
	var ps5_m6;
	var om6_m6;
	var ph6_m6;
	var ps6_m6;
	var om7_m6;
	var ph7_m6;
	var ps7_m6;
	var om8_m6;
	var ph8_m6;
	var ps8_m6;
	var om9_m6;
	var ph9_m6;
	var ps9_m6;
	var om10_m6;
	var ph10_m6;
	var ps10_m6 >= 1.0 ,  <= 1.0;
	var omm10_m5;
	var phm10_m5;
	var psm10_m5 >= 1.0 ,  <= 1.0;
	var omm9_m5;
	var phm9_m5;
	var psm9_m5;
	var omm8_m5;
	var phm8_m5;
	var psm8_m5;
	var omm7_m5;
	var phm7_m5;
	var psm7_m5;
	var omm6_m5;
	var phm6_m5;
	var psm6_m5;
	var omm5_m5;
	var phm5_m5;
	var psm5_m5;
	var omm4_m5;
	var phm4_m5;
	var psm4_m5;
	var omm3_m5;
	var phm3_m5;
	var psm3_m5;
	var omm2_m5;
	var phm2_m5;
	var psm2_m5;
	var omm1_m5;
	var phm1_m5;
	var psm1_m5;
	var om0_m5;
	var ph0_m5;
	var ps0_m5;
	var om1_m5;
	var ph1_m5;
	var ps1_m5;
	var om2_m5;
	var ph2_m5;
	var ps2_m5;
	var om3_m5;
	var ph3_m5;
	var ps3_m5;
	var om4_m5;
	var ph4_m5;
	var ps4_m5;
	var om5_m5;
	var ph5_m5;
	var ps5_m5;
	var om6_m5;
	var ph6_m5;
	var ps6_m5;
	var om7_m5;
	var ph7_m5;
	var ps7_m5;
	var om8_m5;
	var ph8_m5;
	var ps8_m5;
	var om9_m5;
	var ph9_m5;
	var ps9_m5;
	var om10_m5;
	var ph10_m5;
	var ps10_m5 >= 1.0 ,  <= 1.0;
	var omm10_m4;
	var phm10_m4;
	var psm10_m4 >= 1.0 ,  <= 1.0;
	var omm9_m4;
	var phm9_m4;
	var psm9_m4;
	var omm8_m4;
	var phm8_m4;
	var psm8_m4;
	var omm7_m4;
	var phm7_m4;
	var psm7_m4;
	var omm6_m4;
	var phm6_m4;
	var psm6_m4;
	var omm5_m4;
	var phm5_m4;
	var psm5_m4;
	var omm4_m4;
	var phm4_m4;
	var psm4_m4;
	var omm3_m4;
	var phm3_m4;
	var psm3_m4;
	var omm2_m4;
	var phm2_m4;
	var psm2_m4;
	var omm1_m4;
	var phm1_m4;
	var psm1_m4;
	var om0_m4;
	var ph0_m4;
	var ps0_m4;
	var om1_m4;
	var ph1_m4;
	var ps1_m4;
	var om2_m4;
	var ph2_m4;
	var ps2_m4;
	var om3_m4;
	var ph3_m4;
	var ps3_m4;
	var om4_m4;
	var ph4_m4;
	var ps4_m4;
	var om5_m4;
	var ph5_m4;
	var ps5_m4;
	var om6_m4;
	var ph6_m4;
	var ps6_m4;
	var om7_m4;
	var ph7_m4;
	var ps7_m4;
	var om8_m4;
	var ph8_m4;
	var ps8_m4;
	var om9_m4;
	var ph9_m4;
	var ps9_m4;
	var om10_m4;
	var ph10_m4;
	var ps10_m4 >= 1.0 ,  <= 1.0;
	var omm10_m3;
	var phm10_m3;
	var psm10_m3 >= 1.0 ,  <= 1.0;
	var omm9_m3;
	var phm9_m3;
	var psm9_m3;
	var omm8_m3;
	var phm8_m3;
	var psm8_m3;
	var omm7_m3;
	var phm7_m3;
	var psm7_m3;
	var omm6_m3;
	var phm6_m3;
	var psm6_m3;
	var omm5_m3;
	var phm5_m3;
	var psm5_m3;
	var omm4_m3;
	var phm4_m3;
	var psm4_m3;
	var omm3_m3;
	var phm3_m3;
	var psm3_m3;
	var omm2_m3;
	var phm2_m3;
	var psm2_m3;
	var omm1_m3;
	var phm1_m3;
	var psm1_m3;
	var om0_m3;
	var ph0_m3;
	var ps0_m3;
	var om1_m3;
	var ph1_m3;
	var ps1_m3;
	var om2_m3;
	var ph2_m3;
	var ps2_m3;
	var om3_m3;
	var ph3_m3;
	var ps3_m3;
	var om4_m3;
	var ph4_m3;
	var ps4_m3;
	var om5_m3;
	var ph5_m3;
	var ps5_m3;
	var om6_m3;
	var ph6_m3;
	var ps6_m3;
	var om7_m3;
	var ph7_m3;
	var ps7_m3;
	var om8_m3;
	var ph8_m3;
	var ps8_m3;
	var om9_m3;
	var ph9_m3;
	var ps9_m3;
	var om10_m3;
	var ph10_m3;
	var ps10_m3 >= 1.0 ,  <= 1.0;
	var omm10_m2;
	var phm10_m2;
	var psm10_m2 >= 1.0 ,  <= 1.0;
	var omm9_m2;
	var phm9_m2;
	var psm9_m2;
	var omm8_m2;
	var phm8_m2;
	var psm8_m2;
	var omm7_m2;
	var phm7_m2;
	var psm7_m2;
	var omm6_m2;
	var phm6_m2;
	var psm6_m2;
	var omm5_m2;
	var phm5_m2;
	var psm5_m2;
	var omm4_m2;
	var phm4_m2;
	var psm4_m2;
	var omm3_m2;
	var phm3_m2;
	var psm3_m2;
	var omm2_m2;
	var phm2_m2;
	var psm2_m2;
	var omm1_m2;
	var phm1_m2;
	var psm1_m2;
	var om0_m2;
	var ph0_m2;
	var ps0_m2;
	var om1_m2;
	var ph1_m2;
	var ps1_m2;
	var om2_m2;
	var ph2_m2;
	var ps2_m2;
	var om3_m2;
	var ph3_m2;
	var ps3_m2;
	var om4_m2;
	var ph4_m2;
	var ps4_m2;
	var om5_m2;
	var ph5_m2;
	var ps5_m2;
	var om6_m2;
	var ph6_m2;
	var ps6_m2;
	var om7_m2;
	var ph7_m2;
	var ps7_m2;
	var om8_m2;
	var ph8_m2;
	var ps8_m2;
	var om9_m2;
	var ph9_m2;
	var ps9_m2;
	var om10_m2;
	var ph10_m2;
	var ps10_m2 >= 1.0 ,  <= 1.0;
	var omm10_m1;
	var phm10_m1;
	var psm10_m1 >= 1.0 ,  <= 1.0;
	var omm9_m1;
	var phm9_m1;
	var psm9_m1;
	var omm8_m1;
	var phm8_m1;
	var psm8_m1;
	var omm7_m1;
	var phm7_m1;
	var psm7_m1;
	var omm6_m1;
	var phm6_m1;
	var psm6_m1;
	var omm5_m1;
	var phm5_m1;
	var psm5_m1;
	var omm4_m1;
	var phm4_m1;
	var psm4_m1;
	var omm3_m1;
	var phm3_m1;
	var psm3_m1;
	var omm2_m1;
	var phm2_m1;
	var psm2_m1;
	var omm1_m1;
	var phm1_m1;
	var psm1_m1;
	var om0_m1;
	var ph0_m1;
	var ps0_m1;
	var om1_m1;
	var ph1_m1;
	var ps1_m1;
	var om2_m1;
	var ph2_m1;
	var ps2_m1;
	var om3_m1;
	var ph3_m1;
	var ps3_m1;
	var om4_m1;
	var ph4_m1;
	var ps4_m1;
	var om5_m1;
	var ph5_m1;
	var ps5_m1;
	var om6_m1;
	var ph6_m1;
	var ps6_m1;
	var om7_m1;
	var ph7_m1;
	var ps7_m1;
	var om8_m1;
	var ph8_m1;
	var ps8_m1;
	var om9_m1;
	var ph9_m1;
	var ps9_m1;
	var om10_m1;
	var ph10_m1;
	var ps10_m1 >= 1.0 ,  <= 1.0;
	var omm10_0;
	var phm10_0;
	var psm10_0 >= 1.0 ,  <= 1.0;
	var omm9_0;
	var phm9_0;
	var psm9_0;
	var omm8_0;
	var phm8_0;
	var psm8_0;
	var omm7_0;
	var phm7_0;
	var psm7_0;
	var omm6_0;
	var phm6_0;
	var psm6_0;
	var omm5_0;
	var phm5_0;
	var psm5_0;
	var omm4_0;
	var phm4_0;
	var psm4_0;
	var omm3_0;
	var phm3_0;
	var psm3_0;
	var omm2_0;
	var phm2_0;
	var psm2_0;
	var omm1_0;
	var phm1_0;
	var psm1_0;
	var om0_0;
	var ph0_0;
	var ps0_0;
	var om1_0;
	var ph1_0;
	var ps1_0;
	var om2_0;
	var ph2_0;
	var ps2_0;
	var om3_0;
	var ph3_0;
	var ps3_0;
	var om4_0;
	var ph4_0;
	var ps4_0;
	var om5_0;
	var ph5_0;
	var ps5_0;
	var om6_0;
	var ph6_0;
	var ps6_0;
	var om7_0;
	var ph7_0;
	var ps7_0;
	var om8_0;
	var ph8_0;
	var ps8_0;
	var om9_0;
	var ph9_0;
	var ps9_0;
	var om10_0;
	var ph10_0;
	var ps10_0 >= 1.0 ,  <= 1.0;
	var omm10_1;
	var phm10_1;
	var psm10_1 >= 1.0 ,  <= 1.0;
	var omm9_1;
	var phm9_1;
	var psm9_1;
	var omm8_1;
	var phm8_1;
	var psm8_1;
	var omm7_1;
	var phm7_1;
	var psm7_1;
	var omm6_1;
	var phm6_1;
	var psm6_1;
	var omm5_1;
	var phm5_1;
	var psm5_1;
	var omm4_1;
	var phm4_1;
	var psm4_1;
	var omm3_1;
	var phm3_1;
	var psm3_1;
	var omm2_1;
	var phm2_1;
	var psm2_1;
	var omm1_1;
	var phm1_1;
	var psm1_1;
	var om0_1;
	var ph0_1;
	var ps0_1;
	var om1_1;
	var ph1_1;
	var ps1_1;
	var om2_1;
	var ph2_1;
	var ps2_1;
	var om3_1;
	var ph3_1;
	var ps3_1;
	var om4_1;
	var ph4_1;
	var ps4_1;
	var om5_1;
	var ph5_1;
	var ps5_1;
	var om6_1;
	var ph6_1;
	var ps6_1;
	var om7_1;
	var ph7_1;
	var ps7_1;
	var om8_1;
	var ph8_1;
	var ps8_1;
	var om9_1;
	var ph9_1;
	var ps9_1;
	var om10_1;
	var ph10_1;
	var ps10_1 >= 1.0 ,  <= 1.0;
	var omm10_2;
	var phm10_2;
	var psm10_2 >= 1.0 ,  <= 1.0;
	var omm9_2;
	var phm9_2;
	var psm9_2;
	var omm8_2;
	var phm8_2;
	var psm8_2;
	var omm7_2;
	var phm7_2;
	var psm7_2;
	var omm6_2;
	var phm6_2;
	var psm6_2;
	var omm5_2;
	var phm5_2;
	var psm5_2;
	var omm4_2;
	var phm4_2;
	var psm4_2;
	var omm3_2;
	var phm3_2;
	var psm3_2;
	var omm2_2;
	var phm2_2;
	var psm2_2;
	var omm1_2;
	var phm1_2;
	var psm1_2;
	var om0_2;
	var ph0_2;
	var ps0_2;
	var om1_2;
	var ph1_2;
	var ps1_2;
	var om2_2;
	var ph2_2;
	var ps2_2;
	var om3_2;
	var ph3_2;
	var ps3_2;
	var om4_2;
	var ph4_2;
	var ps4_2;
	var om5_2;
	var ph5_2;
	var ps5_2;
	var om6_2;
	var ph6_2;
	var ps6_2;
	var om7_2;
	var ph7_2;
	var ps7_2;
	var om8_2;
	var ph8_2;
	var ps8_2;
	var om9_2;
	var ph9_2;
	var ps9_2;
	var om10_2;
	var ph10_2;
	var ps10_2 >= 1.0 ,  <= 1.0;
	var omm10_3;
	var phm10_3;
	var psm10_3 >= 1.0 ,  <= 1.0;
	var omm9_3;
	var phm9_3;
	var psm9_3;
	var omm8_3;
	var phm8_3;
	var psm8_3;
	var omm7_3;
	var phm7_3;
	var psm7_3;
	var omm6_3;
	var phm6_3;
	var psm6_3;
	var omm5_3;
	var phm5_3;
	var psm5_3;
	var omm4_3;
	var phm4_3;
	var psm4_3;
	var omm3_3;
	var phm3_3;
	var psm3_3;
	var omm2_3;
	var phm2_3;
	var psm2_3;
	var omm1_3;
	var phm1_3;
	var psm1_3;
	var om0_3;
	var ph0_3;
	var ps0_3;
	var om1_3;
	var ph1_3;
	var ps1_3;
	var om2_3;
	var ph2_3;
	var ps2_3;
	var om3_3;
	var ph3_3;
	var ps3_3;
	var om4_3;
	var ph4_3;
	var ps4_3;
	var om5_3;
	var ph5_3;
	var ps5_3;
	var om6_3;
	var ph6_3;
	var ps6_3;
	var om7_3;
	var ph7_3;
	var ps7_3;
	var om8_3;
	var ph8_3;
	var ps8_3;
	var om9_3;
	var ph9_3;
	var ps9_3;
	var om10_3;
	var ph10_3;
	var ps10_3 >= 1.0 ,  <= 1.0;
	var omm10_4;
	var phm10_4;
	var psm10_4 >= 1.0 ,  <= 1.0;
	var omm9_4;
	var phm9_4;
	var psm9_4;
	var omm8_4;
	var phm8_4;
	var psm8_4;
	var omm7_4;
	var phm7_4;
	var psm7_4;
	var omm6_4;
	var phm6_4;
	var psm6_4;
	var omm5_4;
	var phm5_4;
	var psm5_4;
	var omm4_4;
	var phm4_4;
	var psm4_4;
	var omm3_4;
	var phm3_4;
	var psm3_4;
	var omm2_4;
	var phm2_4;
	var psm2_4;
	var omm1_4;
	var phm1_4;
	var psm1_4;
	var om0_4;
	var ph0_4;
	var ps0_4;
	var om1_4;
	var ph1_4;
	var ps1_4;
	var om2_4;
	var ph2_4;
	var ps2_4;
	var om3_4;
	var ph3_4;
	var ps3_4;
	var om4_4;
	var ph4_4;
	var ps4_4;
	var om5_4;
	var ph5_4;
	var ps5_4;
	var om6_4;
	var ph6_4;
	var ps6_4;
	var om7_4;
	var ph7_4;
	var ps7_4;
	var om8_4;
	var ph8_4;
	var ps8_4;
	var om9_4;
	var ph9_4;
	var ps9_4;
	var om10_4;
	var ph10_4;
	var ps10_4 >= 1.0 ,  <= 1.0;
	var omm10_5;
	var phm10_5;
	var psm10_5 >= 1.0 ,  <= 1.0;
	var omm9_5;
	var phm9_5;
	var psm9_5;
	var omm8_5;
	var phm8_5;
	var psm8_5;
	var omm7_5;
	var phm7_5;
	var psm7_5;
	var omm6_5;
	var phm6_5;
	var psm6_5;
	var omm5_5;
	var phm5_5;
	var psm5_5;
	var omm4_5;
	var phm4_5;
	var psm4_5;
	var omm3_5;
	var phm3_5;
	var psm3_5;
	var omm2_5;
	var phm2_5;
	var psm2_5;
	var omm1_5;
	var phm1_5;
	var psm1_5;
	var om0_5;
	var ph0_5;
	var ps0_5;
	var om1_5;
	var ph1_5;
	var ps1_5;
	var om2_5;
	var ph2_5;
	var ps2_5;
	var om3_5;
	var ph3_5;
	var ps3_5;
	var om4_5;
	var ph4_5;
	var ps4_5;
	var om5_5;
	var ph5_5;
	var ps5_5;
	var om6_5;
	var ph6_5;
	var ps6_5;
	var om7_5;
	var ph7_5;
	var ps7_5;
	var om8_5;
	var ph8_5;
	var ps8_5;
	var om9_5;
	var ph9_5;
	var ps9_5;
	var om10_5;
	var ph10_5;
	var ps10_5 >= 1.0 ,  <= 1.0;
	var omm10_6;
	var phm10_6;
	var psm10_6 >= 1.0 ,  <= 1.0;
	var omm9_6;
	var phm9_6;
	var psm9_6;
	var omm8_6;
	var phm8_6;
	var psm8_6;
	var omm7_6;
	var phm7_6;
	var psm7_6;
	var omm6_6;
	var phm6_6;
	var psm6_6;
	var omm5_6;
	var phm5_6;
	var psm5_6;
	var omm4_6;
	var phm4_6;
	var psm4_6;
	var omm3_6;
	var phm3_6;
	var psm3_6;
	var omm2_6;
	var phm2_6;
	var psm2_6;
	var omm1_6;
	var phm1_6;
	var psm1_6;
	var om0_6;
	var ph0_6;
	var ps0_6;
	var om1_6;
	var ph1_6;
	var ps1_6;
	var om2_6;
	var ph2_6;
	var ps2_6;
	var om3_6;
	var ph3_6;
	var ps3_6;
	var om4_6;
	var ph4_6;
	var ps4_6;
	var om5_6;
	var ph5_6;
	var ps5_6;
	var om6_6;
	var ph6_6;
	var ps6_6;
	var om7_6;
	var ph7_6;
	var ps7_6;
	var om8_6;
	var ph8_6;
	var ps8_6;
	var om9_6;
	var ph9_6;
	var ps9_6;
	var om10_6;
	var ph10_6;
	var ps10_6 >= 1.0 ,  <= 1.0;
	var omm10_7;
	var phm10_7;
	var psm10_7 >= 1.0 ,  <= 1.0;
	var omm9_7;
	var phm9_7;
	var psm9_7;
	var omm8_7;
	var phm8_7;
	var psm8_7;
	var omm7_7;
	var phm7_7;
	var psm7_7;
	var omm6_7;
	var phm6_7;
	var psm6_7;
	var omm5_7;
	var phm5_7;
	var psm5_7;
	var omm4_7;
	var phm4_7;
	var psm4_7;
	var omm3_7;
	var phm3_7;
	var psm3_7;
	var omm2_7;
	var phm2_7;
	var psm2_7;
	var omm1_7;
	var phm1_7;
	var psm1_7;
	var om0_7;
	var ph0_7;
	var ps0_7;
	var om1_7;
	var ph1_7;
	var ps1_7;
	var om2_7;
	var ph2_7;
	var ps2_7;
	var om3_7;
	var ph3_7;
	var ps3_7;
	var om4_7;
	var ph4_7;
	var ps4_7;
	var om5_7;
	var ph5_7;
	var ps5_7;
	var om6_7;
	var ph6_7;
	var ps6_7;
	var om7_7;
	var ph7_7;
	var ps7_7;
	var om8_7;
	var ph8_7;
	var ps8_7;
	var om9_7;
	var ph9_7;
	var ps9_7;
	var om10_7;
	var ph10_7;
	var ps10_7 >= 1.0 ,  <= 1.0;
	var omm10_8;
	var phm10_8;
	var psm10_8 >= 1.0 ,  <= 1.0;
	var omm9_8;
	var phm9_8;
	var psm9_8;
	var omm8_8;
	var phm8_8;
	var psm8_8;
	var omm7_8;
	var phm7_8;
	var psm7_8;
	var omm6_8;
	var phm6_8;
	var psm6_8;
	var omm5_8;
	var phm5_8;
	var psm5_8;
	var omm4_8;
	var phm4_8;
	var psm4_8;
	var omm3_8;
	var phm3_8;
	var psm3_8;
	var omm2_8;
	var phm2_8;
	var psm2_8;
	var omm1_8;
	var phm1_8;
	var psm1_8;
	var om0_8;
	var ph0_8;
	var ps0_8;
	var om1_8;
	var ph1_8;
	var ps1_8;
	var om2_8;
	var ph2_8;
	var ps2_8;
	var om3_8;
	var ph3_8;
	var ps3_8;
	var om4_8;
	var ph4_8;
	var ps4_8;
	var om5_8;
	var ph5_8;
	var ps5_8;
	var om6_8;
	var ph6_8;
	var ps6_8;
	var om7_8;
	var ph7_8;
	var ps7_8;
	var om8_8;
	var ph8_8;
	var ps8_8;
	var om9_8;
	var ph9_8;
	var ps9_8;
	var om10_8;
	var ph10_8;
	var ps10_8 >= 1.0 ,  <= 1.0;
	var omm10_9;
	var phm10_9;
	var psm10_9 >= 1.0 ,  <= 1.0;
	var omm9_9;
	var phm9_9;
	var psm9_9;
	var omm8_9;
	var phm8_9;
	var psm8_9;
	var omm7_9;
	var phm7_9;
	var psm7_9;
	var omm6_9;
	var phm6_9;
	var psm6_9;
	var omm5_9;
	var phm5_9;
	var psm5_9;
	var omm4_9;
	var phm4_9;
	var psm4_9;
	var omm3_9;
	var phm3_9;
	var psm3_9;
	var omm2_9;
	var phm2_9;
	var psm2_9;
	var omm1_9;
	var phm1_9;
	var psm1_9;
	var om0_9;
	var ph0_9;
	var ps0_9;
	var om1_9;
	var ph1_9;
	var ps1_9;
	var om2_9;
	var ph2_9;
	var ps2_9;
	var om3_9;
	var ph3_9;
	var ps3_9;
	var om4_9;
	var ph4_9;
	var ps4_9;
	var om5_9;
	var ph5_9;
	var ps5_9;
	var om6_9;
	var ph6_9;
	var ps6_9;
	var om7_9;
	var ph7_9;
	var ps7_9;
	var om8_9;
	var ph8_9;
	var ps8_9;
	var om9_9;
	var ph9_9;
	var ps9_9;
	var om10_9;
	var ph10_9;
	var ps10_9 >= 1.0 ,  <= 1.0;
	var omm10_10;
	var phm10_10;
	var psm10_10 >= 1.0 ,  <= 1.0;
	var omm9_10;
	var phm9_10;
	var psm9_10 >= 1.0 ,  <= 1.0;
	var omm8_10;
	var phm8_10;
	var psm8_10 >= 1.0 ,  <= 1.0;
	var omm7_10;
	var phm7_10;
	var psm7_10 >= 1.0 ,  <= 1.0;
	var omm6_10;
	var phm6_10;
	var psm6_10 >= 1.0 ,  <= 1.0;
	var omm5_10;
	var phm5_10;
	var psm5_10 >= 1.0 ,  <= 1.0;
	var omm4_10;
	var phm4_10;
	var psm4_10 >= 1.0 ,  <= 1.0;
	var omm3_10;
	var phm3_10;
	var psm3_10 >= 1.0 ,  <= 1.0;
	var omm2_10;
	var phm2_10;
	var psm2_10 >= 1.0 ,  <= 1.0;
	var omm1_10;
	var phm1_10;
	var psm1_10 >= 1.0 ,  <= 1.0;
	var om0_10;
	var ph0_10;
	var ps0_10 >= 1.0 ,  <= 1.0;
	var om1_10;
	var ph1_10;
	var ps1_10 >= 1.0 ,  <= 1.0;
	var om2_10;
	var ph2_10;
	var ps2_10 >= 1.0 ,  <= 1.0;
	var om3_10;
	var ph3_10;
	var ps3_10 >= 1.0 ,  <= 1.0;
	var om4_10;
	var ph4_10;
	var ps4_10 >= 1.0 ,  <= 1.0;
	var om5_10;
	var ph5_10;
	var ps5_10 >= 1.0 ,  <= 1.0;
	var om6_10;
	var ph6_10;
	var ps6_10 >= 1.0 ,  <= 1.0;
	var om7_10;
	var ph7_10;
	var ps7_10 >= 1.0 ,  <= 1.0;
	var om8_10;
	var ph8_10;
	var ps8_10 >= 1.0 ,  <= 1.0;
	var om9_10;
	var ph9_10;
	var ps9_10 >= 1.0 ,  <= 1.0;
	var om10_10;
	var ph10_10;
	var ps10_10 >= 1.0 ,  <= 1.0;

minimize obj:
	(-200.0*omm9_9 + 100.0*omm8_9 + 100.0*omm10_9 - 200.0*omm9_9 + 
	100.0*omm9_8 + 100.0*omm9_10 + 1.5308084989341916e-13*phm8_9 - 
	1.5308084989341916e-13*phm10_9 - 2500.0*phm9_8 + 2500.0*phm9_10)^2 + 
	(-200.0*psm9_9 + 100.0*psm8_9 + 100.0*psm10_9 - 200.0*psm9_9 + 
	100.0*psm9_8 + 100.0*psm9_10 + 0.25*omm9_9)^2 + (-25.0*(psm9_8-psm9_10) * 
	(phm8_9-phm10_9) + 25.0*(psm8_9-psm10_9) * (phm9_8-phm9_10) - 
	200.0*phm9_9 + 100.0*phm8_9 + 100.0*phm10_9 - 200.0*phm9_9 + 100.0*phm9_8 
	+ 100.0*phm9_10)^2 + (-200.0*omm8_9 + 100.0*omm7_9 + 100.0*omm9_9 - 
	200.0*omm8_9 + 100.0*omm8_8 + 100.0*omm8_10 + 1.5308084989341916e-13*phm7_9 
	- 1.5308084989341916e-13*phm9_9 - 2500.0*phm8_8 + 2500.0*phm8_10)^2 + 
	(-200.0*psm8_9 + 100.0*psm7_9 + 100.0*psm9_9 - 200.0*psm8_9 + 100.0*psm8_8 
	+ 100.0*psm8_10 + 0.25*omm8_9)^2 + (-25.0*(psm8_8-psm8_10) * 
	(phm7_9-phm9_9) + 25.0*(psm7_9-psm9_9) * (phm8_8-phm8_10) - 200.0*phm8_9 
	+ 100.0*phm7_9 + 100.0*phm9_9 - 200.0*phm8_9 + 100.0*phm8_8 + 
	100.0*phm8_10)^2 + (-200.0*omm7_9 + 100.0*omm6_9 + 100.0*omm8_9 - 
	200.0*omm7_9 + 100.0*omm7_8 + 100.0*omm7_10 + 1.5308084989341916e-13*phm6_9 
	- 1.5308084989341916e-13*phm8_9 - 2500.0*phm7_8 + 2500.0*phm7_10)^2 + 
	(-200.0*psm7_9 + 100.0*psm6_9 + 100.0*psm8_9 - 200.0*psm7_9 + 100.0*psm7_8 
	+ 100.0*psm7_10 + 0.25*omm7_9)^2 + (-25.0*(psm7_8-psm7_10) * 
	(phm6_9-phm8_9) + 25.0*(psm6_9-psm8_9) * (phm7_8-phm7_10) - 200.0*phm7_9 
	+ 100.0*phm6_9 + 100.0*phm8_9 - 200.0*phm7_9 + 100.0*phm7_8 + 
	100.0*phm7_10)^2 + (-200.0*omm6_9 + 100.0*omm5_9 + 100.0*omm7_9 - 
	200.0*omm6_9 + 100.0*omm6_8 + 100.0*omm6_10 + 1.5308084989341916e-13*phm5_9 
	- 1.5308084989341916e-13*phm7_9 - 2500.0*phm6_8 + 2500.0*phm6_10)^2 + 
	(-200.0*psm6_9 + 100.0*psm5_9 + 100.0*psm7_9 - 200.0*psm6_9 + 100.0*psm6_8 
	+ 100.0*psm6_10 + 0.25*omm6_9)^2 + (-25.0*(psm6_8-psm6_10) * 
	(phm5_9-phm7_9) + 25.0*(psm5_9-psm7_9) * (phm6_8-phm6_10) - 200.0*phm6_9 
	+ 100.0*phm5_9 + 100.0*phm7_9 - 200.0*phm6_9 + 100.0*phm6_8 + 
	100.0*phm6_10)^2 + (-200.0*omm5_9 + 100.0*omm4_9 + 100.0*omm6_9 - 
	200.0*omm5_9 + 100.0*omm5_8 + 100.0*omm5_10 + 1.5308084989341916e-13*phm4_9 
	- 1.5308084989341916e-13*phm6_9 - 2500.0*phm5_8 + 2500.0*phm5_10)^2 + 
	(-200.0*psm5_9 + 100.0*psm4_9 + 100.0*psm6_9 - 200.0*psm5_9 + 100.0*psm5_8 
	+ 100.0*psm5_10 + 0.25*omm5_9)^2 + (-25.0*(psm5_8-psm5_10) * 
	(phm4_9-phm6_9) + 25.0*(psm4_9-psm6_9) * (phm5_8-phm5_10) - 200.0*phm5_9 
	+ 100.0*phm4_9 + 100.0*phm6_9 - 200.0*phm5_9 + 100.0*phm5_8 + 
	100.0*phm5_10)^2 + (-200.0*omm4_9 + 100.0*omm3_9 + 100.0*omm5_9 - 
	200.0*omm4_9 + 100.0*omm4_8 + 100.0*omm4_10 + 1.5308084989341916e-13*phm3_9 
	- 1.5308084989341916e-13*phm5_9 - 2500.0*phm4_8 + 2500.0*phm4_10)^2 + 
	(-200.0*psm4_9 + 100.0*psm3_9 + 100.0*psm5_9 - 200.0*psm4_9 + 100.0*psm4_8 
	+ 100.0*psm4_10 + 0.25*omm4_9)^2 + (-25.0*(psm4_8-psm4_10) * 
	(phm3_9-phm5_9) + 25.0*(psm3_9-psm5_9) * (phm4_8-phm4_10) - 200.0*phm4_9 
	+ 100.0*phm3_9 + 100.0*phm5_9 - 200.0*phm4_9 + 100.0*phm4_8 + 
	100.0*phm4_10)^2 + (-200.0*omm3_9 + 100.0*omm2_9 + 100.0*omm4_9 - 
	200.0*omm3_9 + 100.0*omm3_8 + 100.0*omm3_10 + 1.5308084989341916e-13*phm2_9 
	- 1.5308084989341916e-13*phm4_9 - 2500.0*phm3_8 + 2500.0*phm3_10)^2 + 
	(-200.0*psm3_9 + 100.0*psm2_9 + 100.0*psm4_9 - 200.0*psm3_9 + 100.0*psm3_8 
	+ 100.0*psm3_10 + 0.25*omm3_9)^2 + (-25.0*(psm3_8-psm3_10) * 
	(phm2_9-phm4_9) + 25.0*(psm2_9-psm4_9) * (phm3_8-phm3_10) - 200.0*phm3_9 
	+ 100.0*phm2_9 + 100.0*phm4_9 - 200.0*phm3_9 + 100.0*phm3_8 + 
	100.0*phm3_10)^2 + (-200.0*omm2_9 + 100.0*omm1_9 + 100.0*omm3_9 - 
	200.0*omm2_9 + 100.0*omm2_8 + 100.0*omm2_10 + 1.5308084989341916e-13*phm1_9 
	- 1.5308084989341916e-13*phm3_9 - 2500.0*phm2_8 + 2500.0*phm2_10)^2 + 
	(-200.0*psm2_9 + 100.0*psm1_9 + 100.0*psm3_9 - 200.0*psm2_9 + 100.0*psm2_8 
	+ 100.0*psm2_10 + 0.25*omm2_9)^2 + (-25.0*(psm2_8-psm2_10) * 
	(phm1_9-phm3_9) + 25.0*(psm1_9-psm3_9) * (phm2_8-phm2_10) - 200.0*phm2_9 
	+ 100.0*phm1_9 + 100.0*phm3_9 - 200.0*phm2_9 + 100.0*phm2_8 + 
	100.0*phm2_10)^2 + (-200.0*omm1_9 + 100.0*om0_9 + 100.0*omm2_9 - 
	200.0*omm1_9 + 100.0*omm1_8 + 100.0*omm1_10 + 1.5308084989341916e-13*ph0_9 
	- 1.5308084989341916e-13*phm2_9 - 2500.0*phm1_8 + 2500.0*phm1_10)^2 + 
	(-200.0*psm1_9 + 100.0*ps0_9 + 100.0*psm2_9 - 200.0*psm1_9 + 100.0*psm1_8 
	+ 100.0*psm1_10 + 0.25*omm1_9)^2 + (-25.0*(psm1_8-psm1_10) * 
	(ph0_9-phm2_9) + 25.0*(ps0_9-psm2_9) * (phm1_8-phm1_10) - 200.0*phm1_9 + 
	100.0*ph0_9 + 100.0*phm2_9 - 200.0*phm1_9 + 100.0*phm1_8 + 
	100.0*phm1_10)^2 + (-200.0*om0_9 + 100.0*om1_9 + 100.0*omm1_9 - 
	200.0*om0_9 + 100.0*om0_8 + 100.0*om0_10 + 1.5308084989341916e-13*ph1_9 - 
	1.5308084989341916e-13*phm1_9 - 2500.0*ph0_8 + 2500.0*ph0_10)^2 + 
	(-200.0*ps0_9 + 100.0*ps1_9 + 100.0*psm1_9 - 200.0*ps0_9 + 100.0*ps0_8 + 
	100.0*ps0_10 + 0.25*om0_9)^2 + (-25.0*(ps0_8-ps0_10) * (ph1_9-phm1_9) + 
	25.0*(ps1_9-psm1_9) * (ph0_8-ph0_10) - 200.0*ph0_9 + 100.0*ph1_9 + 
	100.0*phm1_9 - 200.0*ph0_9 + 100.0*ph0_8 + 100.0*ph0_10)^2 + (-200.0*om1_9 
	+ 100.0*om2_9 + 100.0*om0_9 - 200.0*om1_9 + 100.0*om1_8 + 100.0*om1_10 + 
	1.5308084989341916e-13*ph2_9 - 1.5308084989341916e-13*ph0_9 - 2500.0*ph1_8 + 
	2500.0*ph1_10)^2 + (-200.0*ps1_9 + 100.0*ps2_9 + 100.0*ps0_9 - 200.0*ps1_9 
	+ 100.0*ps1_8 + 100.0*ps1_10 + 0.25*om1_9)^2 + (-25.0*(ps1_8-ps1_10) * 
	(ph2_9-ph0_9) + 25.0*(ps2_9-ps0_9) * (ph1_8-ph1_10) - 200.0*ph1_9 + 
	100.0*ph2_9 + 100.0*ph0_9 - 200.0*ph1_9 + 100.0*ph1_8 + 100.0*ph1_10)^2 + 
	(-200.0*om2_9 + 100.0*om3_9 + 100.0*om1_9 - 200.0*om2_9 + 100.0*om2_8 + 
	100.0*om2_10 + 1.5308084989341916e-13*ph3_9 - 1.5308084989341916e-13*ph1_9 - 
	2500.0*ph2_8 + 2500.0*ph2_10)^2 + (-200.0*ps2_9 + 100.0*ps3_9 + 
	100.0*ps1_9 - 200.0*ps2_9 + 100.0*ps2_8 + 100.0*ps2_10 + 0.25*om2_9)^2 + 
	(-25.0*(ps2_8-ps2_10) * (ph3_9-ph1_9) + 25.0*(ps3_9-ps1_9) * 
	(ph2_8-ph2_10) - 200.0*ph2_9 + 100.0*ph3_9 + 100.0*ph1_9 - 200.0*ph2_9 + 
	100.0*ph2_8 + 100.0*ph2_10)^2 + (-200.0*om3_9 + 100.0*om4_9 + 100.0*om2_9 
	- 200.0*om3_9 + 100.0*om3_8 + 100.0*om3_10 + 1.5308084989341916e-13*ph4_9 - 
	1.5308084989341916e-13*ph2_9 - 2500.0*ph3_8 + 2500.0*ph3_10)^2 + 
	(-200.0*ps3_9 + 100.0*ps4_9 + 100.0*ps2_9 - 200.0*ps3_9 + 100.0*ps3_8 + 
	100.0*ps3_10 + 0.25*om3_9)^2 + (-25.0*(ps3_8-ps3_10) * (ph4_9-ph2_9) + 
	25.0*(ps4_9-ps2_9) * (ph3_8-ph3_10) - 200.0*ph3_9 + 100.0*ph4_9 + 
	100.0*ph2_9 - 200.0*ph3_9 + 100.0*ph3_8 + 100.0*ph3_10)^2 + (-200.0*om4_9 
	+ 100.0*om5_9 + 100.0*om3_9 - 200.0*om4_9 + 100.0*om4_8 + 100.0*om4_10 + 
	1.5308084989341916e-13*ph5_9 - 1.5308084989341916e-13*ph3_9 - 2500.0*ph4_8 + 
	2500.0*ph4_10)^2 + (-200.0*ps4_9 + 100.0*ps5_9 + 100.0*ps3_9 - 200.0*ps4_9 
	+ 100.0*ps4_8 + 100.0*ps4_10 + 0.25*om4_9)^2 + (-25.0*(ps4_8-ps4_10) * 
	(ph5_9-ph3_9) + 25.0*(ps5_9-ps3_9) * (ph4_8-ph4_10) - 200.0*ph4_9 + 
	100.0*ph5_9 + 100.0*ph3_9 - 200.0*ph4_9 + 100.0*ph4_8 + 100.0*ph4_10)^2 + 
	(-200.0*om5_9 + 100.0*om6_9 + 100.0*om4_9 - 200.0*om5_9 + 100.0*om5_8 + 
	100.0*om5_10 + 1.5308084989341916e-13*ph6_9 - 1.5308084989341916e-13*ph4_9 - 
	2500.0*ph5_8 + 2500.0*ph5_10)^2 + (-200.0*ps5_9 + 100.0*ps6_9 + 
	100.0*ps4_9 - 200.0*ps5_9 + 100.0*ps5_8 + 100.0*ps5_10 + 0.25*om5_9)^2 + 
	(-25.0*(ps5_8-ps5_10) * (ph6_9-ph4_9) + 25.0*(ps6_9-ps4_9) * 
	(ph5_8-ph5_10) - 200.0*ph5_9 + 100.0*ph6_9 + 100.0*ph4_9 - 200.0*ph5_9 + 
	100.0*ph5_8 + 100.0*ph5_10)^2 + (-200.0*om6_9 + 100.0*om7_9 + 100.0*om5_9 
	- 200.0*om6_9 + 100.0*om6_8 + 100.0*om6_10 + 1.5308084989341916e-13*ph7_9 - 
	1.5308084989341916e-13*ph5_9 - 2500.0*ph6_8 + 2500.0*ph6_10)^2 + 
	(-200.0*ps6_9 + 100.0*ps7_9 + 100.0*ps5_9 - 200.0*ps6_9 + 100.0*ps6_8 + 
	100.0*ps6_10 + 0.25*om6_9)^2 + (-25.0*(ps6_8-ps6_10) * (ph7_9-ph5_9) + 
	25.0*(ps7_9-ps5_9) * (ph6_8-ph6_10) - 200.0*ph6_9 + 100.0*ph7_9 + 
	100.0*ph5_9 - 200.0*ph6_9 + 100.0*ph6_8 + 100.0*ph6_10)^2 + (-200.0*om7_9 
	+ 100.0*om8_9 + 100.0*om6_9 - 200.0*om7_9 + 100.0*om7_8 + 100.0*om7_10 + 
	1.5308084989341916e-13*ph8_9 - 1.5308084989341916e-13*ph6_9 - 2500.0*ph7_8 + 
	2500.0*ph7_10)^2 + (-200.0*ps7_9 + 100.0*ps8_9 + 100.0*ps6_9 - 200.0*ps7_9 
	+ 100.0*ps7_8 + 100.0*ps7_10 + 0.25*om7_9)^2 + (-25.0*(ps7_8-ps7_10) * 
	(ph8_9-ph6_9) + 25.0*(ps8_9-ps6_9) * (ph7_8-ph7_10) - 200.0*ph7_9 + 
	100.0*ph8_9 + 100.0*ph6_9 - 200.0*ph7_9 + 100.0*ph7_8 + 100.0*ph7_10)^2 + 
	(-200.0*om8_9 + 100.0*om9_9 + 100.0*om7_9 - 200.0*om8_9 + 100.0*om8_8 + 
	100.0*om8_10 + 1.5308084989341916e-13*ph9_9 - 1.5308084989341916e-13*ph7_9 - 
	2500.0*ph8_8 + 2500.0*ph8_10)^2 + (-200.0*ps8_9 + 100.0*ps9_9 + 
	100.0*ps7_9 - 200.0*ps8_9 + 100.0*ps8_8 + 100.0*ps8_10 + 0.25*om8_9)^2 + 
	(-25.0*(ps8_8-ps8_10) * (ph9_9-ph7_9) + 25.0*(ps9_9-ps7_9) * 
	(ph8_8-ph8_10) - 200.0*ph8_9 + 100.0*ph9_9 + 100.0*ph7_9 - 200.0*ph8_9 + 
	100.0*ph8_8 + 100.0*ph8_10)^2 + (-200.0*om9_9 + 100.0*om10_9 + 100.0*om8_9 
	- 200.0*om9_9 + 100.0*om9_8 + 100.0*om9_10 + 1.5308084989341916e-13*ph10_9 
	- 1.5308084989341916e-13*ph8_9 - 2500.0*ph9_8 + 2500.0*ph9_10)^2 + 
	(-200.0*ps9_9 + 100.0*ps10_9 + 100.0*ps8_9 - 200.0*ps9_9 + 100.0*ps9_8 + 
	100.0*ps9_10 + 0.25*om9_9)^2 + (-25.0*(ps9_8-ps9_10) * (ph10_9-ph8_9) + 
	25.0*(ps10_9-ps8_9) * (ph9_8-ph9_10) - 200.0*ph9_9 + 100.0*ph10_9 + 
	100.0*ph8_9 - 200.0*ph9_9 + 100.0*ph9_8 + 100.0*ph9_10)^2 + (-200.0*omm9_8 
	+ 100.0*omm8_8 + 100.0*omm10_8 - 200.0*omm9_8 + 100.0*omm9_7 + 
	100.0*omm9_9 + 1.5308084989341916e-13*phm8_8 - 
	1.5308084989341916e-13*phm10_8 - 2500.0*phm9_7 + 2500.0*phm9_9)^2 + 
	(-200.0*psm9_8 + 100.0*psm8_8 + 100.0*psm10_8 - 200.0*psm9_8 + 
	100.0*psm9_7 + 100.0*psm9_9 + 0.25*omm9_8)^2 + (-25.0*(psm9_7-psm9_9) * 
	(phm8_8-phm10_8) + 25.0*(psm8_8-psm10_8) * (phm9_7-phm9_9) - 
	200.0*phm9_8 + 100.0*phm8_8 + 100.0*phm10_8 - 200.0*phm9_8 + 100.0*phm9_7 
	+ 100.0*phm9_9)^2 + (-200.0*omm8_8 + 100.0*omm7_8 + 100.0*omm9_8 - 
	200.0*omm8_8 + 100.0*omm8_7 + 100.0*omm8_9 + 1.5308084989341916e-13*phm7_8 
	- 1.5308084989341916e-13*phm9_8 - 2500.0*phm8_7 + 2500.0*phm8_9)^2 + 
	(-200.0*psm8_8 + 100.0*psm7_8 + 100.0*psm9_8 - 200.0*psm8_8 + 100.0*psm8_7 
	+ 100.0*psm8_9 + 0.25*omm8_8)^2 + (-25.0*(psm8_7-psm8_9) * 
	(phm7_8-phm9_8) + 25.0*(psm7_8-psm9_8) * (phm8_7-phm8_9) - 200.0*phm8_8 
	+ 100.0*phm7_8 + 100.0*phm9_8 - 200.0*phm8_8 + 100.0*phm8_7 + 
	100.0*phm8_9)^2 + (-200.0*omm7_8 + 100.0*omm6_8 + 100.0*omm8_8 - 
	200.0*omm7_8 + 100.0*omm7_7 + 100.0*omm7_9 + 1.5308084989341916e-13*phm6_8 
	- 1.5308084989341916e-13*phm8_8 - 2500.0*phm7_7 + 2500.0*phm7_9)^2 + 
	(-200.0*psm7_8 + 100.0*psm6_8 + 100.0*psm8_8 - 200.0*psm7_8 + 100.0*psm7_7 
	+ 100.0*psm7_9 + 0.25*omm7_8)^2 + (-25.0*(psm7_7-psm7_9) * 
	(phm6_8-phm8_8) + 25.0*(psm6_8-psm8_8) * (phm7_7-phm7_9) - 200.0*phm7_8 
	+ 100.0*phm6_8 + 100.0*phm8_8 - 200.0*phm7_8 + 100.0*phm7_7 + 
	100.0*phm7_9)^2 + (-200.0*omm6_8 + 100.0*omm5_8 + 100.0*omm7_8 - 
	200.0*omm6_8 + 100.0*omm6_7 + 100.0*omm6_9 + 1.5308084989341916e-13*phm5_8 
	- 1.5308084989341916e-13*phm7_8 - 2500.0*phm6_7 + 2500.0*phm6_9)^2 + 
	(-200.0*psm6_8 + 100.0*psm5_8 + 100.0*psm7_8 - 200.0*psm6_8 + 100.0*psm6_7 
	+ 100.0*psm6_9 + 0.25*omm6_8)^2 + (-25.0*(psm6_7-psm6_9) * 
	(phm5_8-phm7_8) + 25.0*(psm5_8-psm7_8) * (phm6_7-phm6_9) - 200.0*phm6_8 
	+ 100.0*phm5_8 + 100.0*phm7_8 - 200.0*phm6_8 + 100.0*phm6_7 + 
	100.0*phm6_9)^2 + (-200.0*omm5_8 + 100.0*omm4_8 + 100.0*omm6_8 - 
	200.0*omm5_8 + 100.0*omm5_7 + 100.0*omm5_9 + 1.5308084989341916e-13*phm4_8 
	- 1.5308084989341916e-13*phm6_8 - 2500.0*phm5_7 + 2500.0*phm5_9)^2 + 
	(-200.0*psm5_8 + 100.0*psm4_8 + 100.0*psm6_8 - 200.0*psm5_8 + 100.0*psm5_7 
	+ 100.0*psm5_9 + 0.25*omm5_8)^2 + (-25.0*(psm5_7-psm5_9) * 
	(phm4_8-phm6_8) + 25.0*(psm4_8-psm6_8) * (phm5_7-phm5_9) - 200.0*phm5_8 
	+ 100.0*phm4_8 + 100.0*phm6_8 - 200.0*phm5_8 + 100.0*phm5_7 + 
	100.0*phm5_9)^2 + (-200.0*omm4_8 + 100.0*omm3_8 + 100.0*omm5_8 - 
	200.0*omm4_8 + 100.0*omm4_7 + 100.0*omm4_9 + 1.5308084989341916e-13*phm3_8 
	- 1.5308084989341916e-13*phm5_8 - 2500.0*phm4_7 + 2500.0*phm4_9)^2 + 
	(-200.0*psm4_8 + 100.0*psm3_8 + 100.0*psm5_8 - 200.0*psm4_8 + 100.0*psm4_7 
	+ 100.0*psm4_9 + 0.25*omm4_8)^2 + (-25.0*(psm4_7-psm4_9) * 
	(phm3_8-phm5_8) + 25.0*(psm3_8-psm5_8) * (phm4_7-phm4_9) - 200.0*phm4_8 
	+ 100.0*phm3_8 + 100.0*phm5_8 - 200.0*phm4_8 + 100.0*phm4_7 + 
	100.0*phm4_9)^2 + (-200.0*omm3_8 + 100.0*omm2_8 + 100.0*omm4_8 - 
	200.0*omm3_8 + 100.0*omm3_7 + 100.0*omm3_9 + 1.5308084989341916e-13*phm2_8 
	- 1.5308084989341916e-13*phm4_8 - 2500.0*phm3_7 + 2500.0*phm3_9)^2 + 
	(-200.0*psm3_8 + 100.0*psm2_8 + 100.0*psm4_8 - 200.0*psm3_8 + 100.0*psm3_7 
	+ 100.0*psm3_9 + 0.25*omm3_8)^2 + (-25.0*(psm3_7-psm3_9) * 
	(phm2_8-phm4_8) + 25.0*(psm2_8-psm4_8) * (phm3_7-phm3_9) - 200.0*phm3_8 
	+ 100.0*phm2_8 + 100.0*phm4_8 - 200.0*phm3_8 + 100.0*phm3_7 + 
	100.0*phm3_9)^2 + (-200.0*omm2_8 + 100.0*omm1_8 + 100.0*omm3_8 - 
	200.0*omm2_8 + 100.0*omm2_7 + 100.0*omm2_9 + 1.5308084989341916e-13*phm1_8 
	- 1.5308084989341916e-13*phm3_8 - 2500.0*phm2_7 + 2500.0*phm2_9)^2 + 
	(-200.0*psm2_8 + 100.0*psm1_8 + 100.0*psm3_8 - 200.0*psm2_8 + 100.0*psm2_7 
	+ 100.0*psm2_9 + 0.25*omm2_8)^2 + (-25.0*(psm2_7-psm2_9) * 
	(phm1_8-phm3_8) + 25.0*(psm1_8-psm3_8) * (phm2_7-phm2_9) - 200.0*phm2_8 
	+ 100.0*phm1_8 + 100.0*phm3_8 - 200.0*phm2_8 + 100.0*phm2_7 + 
	100.0*phm2_9)^2 + (-200.0*omm1_8 + 100.0*om0_8 + 100.0*omm2_8 - 
	200.0*omm1_8 + 100.0*omm1_7 + 100.0*omm1_9 + 1.5308084989341916e-13*ph0_8 - 
	1.5308084989341916e-13*phm2_8 - 2500.0*phm1_7 + 2500.0*phm1_9)^2 + 
	(-200.0*psm1_8 + 100.0*ps0_8 + 100.0*psm2_8 - 200.0*psm1_8 + 100.0*psm1_7 
	+ 100.0*psm1_9 + 0.25*omm1_8)^2 + (-25.0*(psm1_7-psm1_9) * (ph0_8-phm2_8) 
	+ 25.0*(ps0_8-psm2_8) * (phm1_7-phm1_9) - 200.0*phm1_8 + 100.0*ph0_8 + 
	100.0*phm2_8 - 200.0*phm1_8 + 100.0*phm1_7 + 100.0*phm1_9)^2 + 
	(-200.0*om0_8 + 100.0*om1_8 + 100.0*omm1_8 - 200.0*om0_8 + 100.0*om0_7 + 
	100.0*om0_9 + 1.5308084989341916e-13*ph1_8 - 1.5308084989341916e-13*phm1_8 - 
	2500.0*ph0_7 + 2500.0*ph0_9)^2 + (-200.0*ps0_8 + 100.0*ps1_8 + 
	100.0*psm1_8 - 200.0*ps0_8 + 100.0*ps0_7 + 100.0*ps0_9 + 0.25*om0_8)^2 + 
	(-25.0*(ps0_7-ps0_9) * (ph1_8-phm1_8) + 25.0*(ps1_8-psm1_8) * 
	(ph0_7-ph0_9) - 200.0*ph0_8 + 100.0*ph1_8 + 100.0*phm1_8 - 200.0*ph0_8 + 
	100.0*ph0_7 + 100.0*ph0_9)^2 + (-200.0*om1_8 + 100.0*om2_8 + 100.0*om0_8 - 
	200.0*om1_8 + 100.0*om1_7 + 100.0*om1_9 + 1.5308084989341916e-13*ph2_8 - 
	1.5308084989341916e-13*ph0_8 - 2500.0*ph1_7 + 2500.0*ph1_9)^2 + 
	(-200.0*ps1_8 + 100.0*ps2_8 + 100.0*ps0_8 - 200.0*ps1_8 + 100.0*ps1_7 + 
	100.0*ps1_9 + 0.25*om1_8)^2 + (-25.0*(ps1_7-ps1_9) * (ph2_8-ph0_8) + 
	25.0*(ps2_8-ps0_8) * (ph1_7-ph1_9) - 200.0*ph1_8 + 100.0*ph2_8 + 
	100.0*ph0_8 - 200.0*ph1_8 + 100.0*ph1_7 + 100.0*ph1_9)^2 + (-200.0*om2_8 + 
	100.0*om3_8 + 100.0*om1_8 - 200.0*om2_8 + 100.0*om2_7 + 100.0*om2_9 + 
	1.5308084989341916e-13*ph3_8 - 1.5308084989341916e-13*ph1_8 - 2500.0*ph2_7 + 
	2500.0*ph2_9)^2 + (-200.0*ps2_8 + 100.0*ps3_8 + 100.0*ps1_8 - 200.0*ps2_8 
	+ 100.0*ps2_7 + 100.0*ps2_9 + 0.25*om2_8)^2 + (-25.0*(ps2_7-ps2_9) * 
	(ph3_8-ph1_8) + 25.0*(ps3_8-ps1_8) * (ph2_7-ph2_9) - 200.0*ph2_8 + 
	100.0*ph3_8 + 100.0*ph1_8 - 200.0*ph2_8 + 100.0*ph2_7 + 100.0*ph2_9)^2 + 
	(-200.0*om3_8 + 100.0*om4_8 + 100.0*om2_8 - 200.0*om3_8 + 100.0*om3_7 + 
	100.0*om3_9 + 1.5308084989341916e-13*ph4_8 - 1.5308084989341916e-13*ph2_8 - 
	2500.0*ph3_7 + 2500.0*ph3_9)^2 + (-200.0*ps3_8 + 100.0*ps4_8 + 100.0*ps2_8 
	- 200.0*ps3_8 + 100.0*ps3_7 + 100.0*ps3_9 + 0.25*om3_8)^2 + 
	(-25.0*(ps3_7-ps3_9) * (ph4_8-ph2_8) + 25.0*(ps4_8-ps2_8) * 
	(ph3_7-ph3_9) - 200.0*ph3_8 + 100.0*ph4_8 + 100.0*ph2_8 - 200.0*ph3_8 + 
	100.0*ph3_7 + 100.0*ph3_9)^2 + (-200.0*om4_8 + 100.0*om5_8 + 100.0*om3_8 - 
	200.0*om4_8 + 100.0*om4_7 + 100.0*om4_9 + 1.5308084989341916e-13*ph5_8 - 
	1.5308084989341916e-13*ph3_8 - 2500.0*ph4_7 + 2500.0*ph4_9)^2 + 
	(-200.0*ps4_8 + 100.0*ps5_8 + 100.0*ps3_8 - 200.0*ps4_8 + 100.0*ps4_7 + 
	100.0*ps4_9 + 0.25*om4_8)^2 + (-25.0*(ps4_7-ps4_9) * (ph5_8-ph3_8) + 
	25.0*(ps5_8-ps3_8) * (ph4_7-ph4_9) - 200.0*ph4_8 + 100.0*ph5_8 + 
	100.0*ph3_8 - 200.0*ph4_8 + 100.0*ph4_7 + 100.0*ph4_9)^2 + (-200.0*om5_8 + 
	100.0*om6_8 + 100.0*om4_8 - 200.0*om5_8 + 100.0*om5_7 + 100.0*om5_9 + 
	1.5308084989341916e-13*ph6_8 - 1.5308084989341916e-13*ph4_8 - 2500.0*ph5_7 + 
	2500.0*ph5_9)^2 + (-200.0*ps5_8 + 100.0*ps6_8 + 100.0*ps4_8 - 200.0*ps5_8 
	+ 100.0*ps5_7 + 100.0*ps5_9 + 0.25*om5_8)^2 + (-25.0*(ps5_7-ps5_9) * 
	(ph6_8-ph4_8) + 25.0*(ps6_8-ps4_8) * (ph5_7-ph5_9) - 200.0*ph5_8 + 
	100.0*ph6_8 + 100.0*ph4_8 - 200.0*ph5_8 + 100.0*ph5_7 + 100.0*ph5_9)^2 + 
	(-200.0*om6_8 + 100.0*om7_8 + 100.0*om5_8 - 200.0*om6_8 + 100.0*om6_7 + 
	100.0*om6_9 + 1.5308084989341916e-13*ph7_8 - 1.5308084989341916e-13*ph5_8 - 
	2500.0*ph6_7 + 2500.0*ph6_9)^2 + (-200.0*ps6_8 + 100.0*ps7_8 + 100.0*ps5_8 
	- 200.0*ps6_8 + 100.0*ps6_7 + 100.0*ps6_9 + 0.25*om6_8)^2 + 
	(-25.0*(ps6_7-ps6_9) * (ph7_8-ph5_8) + 25.0*(ps7_8-ps5_8) * 
	(ph6_7-ph6_9) - 200.0*ph6_8 + 100.0*ph7_8 + 100.0*ph5_8 - 200.0*ph6_8 + 
	100.0*ph6_7 + 100.0*ph6_9)^2 + (-200.0*om7_8 + 100.0*om8_8 + 100.0*om6_8 - 
	200.0*om7_8 + 100.0*om7_7 + 100.0*om7_9 + 1.5308084989341916e-13*ph8_8 - 
	1.5308084989341916e-13*ph6_8 - 2500.0*ph7_7 + 2500.0*ph7_9)^2 + 
	(-200.0*ps7_8 + 100.0*ps8_8 + 100.0*ps6_8 - 200.0*ps7_8 + 100.0*ps7_7 + 
	100.0*ps7_9 + 0.25*om7_8)^2 + (-25.0*(ps7_7-ps7_9) * (ph8_8-ph6_8) + 
	25.0*(ps8_8-ps6_8) * (ph7_7-ph7_9) - 200.0*ph7_8 + 100.0*ph8_8 + 
	100.0*ph6_8 - 200.0*ph7_8 + 100.0*ph7_7 + 100.0*ph7_9)^2 + (-200.0*om8_8 + 
	100.0*om9_8 + 100.0*om7_8 - 200.0*om8_8 + 100.0*om8_7 + 100.0*om8_9 + 
	1.5308084989341916e-13*ph9_8 - 1.5308084989341916e-13*ph7_8 - 2500.0*ph8_7 + 
	2500.0*ph8_9)^2 + (-200.0*ps8_8 + 100.0*ps9_8 + 100.0*ps7_8 - 200.0*ps8_8 
	+ 100.0*ps8_7 + 100.0*ps8_9 + 0.25*om8_8)^2 + (-25.0*(ps8_7-ps8_9) * 
	(ph9_8-ph7_8) + 25.0*(ps9_8-ps7_8) * (ph8_7-ph8_9) - 200.0*ph8_8 + 
	100.0*ph9_8 + 100.0*ph7_8 - 200.0*ph8_8 + 100.0*ph8_7 + 100.0*ph8_9)^2 + 
	(-200.0*om9_8 + 100.0*om10_8 + 100.0*om8_8 - 200.0*om9_8 + 100.0*om9_7 + 
	100.0*om9_9 + 1.5308084989341916e-13*ph10_8 - 1.5308084989341916e-13*ph8_8 - 
	2500.0*ph9_7 + 2500.0*ph9_9)^2 + (-200.0*ps9_8 + 100.0*ps10_8 + 
	100.0*ps8_8 - 200.0*ps9_8 + 100.0*ps9_7 + 100.0*ps9_9 + 0.25*om9_8)^2 + 
	(-25.0*(ps9_7-ps9_9) * (ph10_8-ph8_8) + 25.0*(ps10_8-ps8_8) * 
	(ph9_7-ph9_9) - 200.0*ph9_8 + 100.0*ph10_8 + 100.0*ph8_8 - 200.0*ph9_8 + 
	100.0*ph9_7 + 100.0*ph9_9)^2 + (-200.0*omm9_7 + 100.0*omm8_7 + 
	100.0*omm10_7 - 200.0*omm9_7 + 100.0*omm9_6 + 100.0*omm9_8 + 
	1.5308084989341916e-13*phm8_7 - 1.5308084989341916e-13*phm10_7 - 
	2500.0*phm9_6 + 2500.0*phm9_8)^2 + (-200.0*psm9_7 + 100.0*psm8_7 + 
	100.0*psm10_7 - 200.0*psm9_7 + 100.0*psm9_6 + 100.0*psm9_8 + 
	0.25*omm9_7)^2 + (-25.0*(psm9_6-psm9_8) * (phm8_7-phm10_7) + 
	25.0*(psm8_7-psm10_7) * (phm9_6-phm9_8) - 200.0*phm9_7 + 100.0*phm8_7 + 
	100.0*phm10_7 - 200.0*phm9_7 + 100.0*phm9_6 + 100.0*phm9_8)^2 + 
	(-200.0*omm8_7 + 100.0*omm7_7 + 100.0*omm9_7 - 200.0*omm8_7 + 100.0*omm8_6 
	+ 100.0*omm8_8 + 1.5308084989341916e-13*phm7_7 - 
	1.5308084989341916e-13*phm9_7 - 2500.0*phm8_6 + 2500.0*phm8_8)^2 + 
	(-200.0*psm8_7 + 100.0*psm7_7 + 100.0*psm9_7 - 200.0*psm8_7 + 100.0*psm8_6 
	+ 100.0*psm8_8 + 0.25*omm8_7)^2 + (-25.0*(psm8_6-psm8_8) * 
	(phm7_7-phm9_7) + 25.0*(psm7_7-psm9_7) * (phm8_6-phm8_8) - 200.0*phm8_7 
	+ 100.0*phm7_7 + 100.0*phm9_7 - 200.0*phm8_7 + 100.0*phm8_6 + 
	100.0*phm8_8)^2 + (-200.0*omm7_7 + 100.0*omm6_7 + 100.0*omm8_7 - 
	200.0*omm7_7 + 100.0*omm7_6 + 100.0*omm7_8 + 1.5308084989341916e-13*phm6_7 
	- 1.5308084989341916e-13*phm8_7 - 2500.0*phm7_6 + 2500.0*phm7_8)^2 + 
	(-200.0*psm7_7 + 100.0*psm6_7 + 100.0*psm8_7 - 200.0*psm7_7 + 100.0*psm7_6 
	+ 100.0*psm7_8 + 0.25*omm7_7)^2 + (-25.0*(psm7_6-psm7_8) * 
	(phm6_7-phm8_7) + 25.0*(psm6_7-psm8_7) * (phm7_6-phm7_8) - 200.0*phm7_7 
	+ 100.0*phm6_7 + 100.0*phm8_7 - 200.0*phm7_7 + 100.0*phm7_6 + 
	100.0*phm7_8)^2 + (-200.0*omm6_7 + 100.0*omm5_7 + 100.0*omm7_7 - 
	200.0*omm6_7 + 100.0*omm6_6 + 100.0*omm6_8 + 1.5308084989341916e-13*phm5_7 
	- 1.5308084989341916e-13*phm7_7 - 2500.0*phm6_6 + 2500.0*phm6_8)^2 + 
	(-200.0*psm6_7 + 100.0*psm5_7 + 100.0*psm7_7 - 200.0*psm6_7 + 100.0*psm6_6 
	+ 100.0*psm6_8 + 0.25*omm6_7)^2 + (-25.0*(psm6_6-psm6_8) * 
	(phm5_7-phm7_7) + 25.0*(psm5_7-psm7_7) * (phm6_6-phm6_8) - 200.0*phm6_7 
	+ 100.0*phm5_7 + 100.0*phm7_7 - 200.0*phm6_7 + 100.0*phm6_6 + 
	100.0*phm6_8)^2 + (-200.0*omm5_7 + 100.0*omm4_7 + 100.0*omm6_7 - 
	200.0*omm5_7 + 100.0*omm5_6 + 100.0*omm5_8 + 1.5308084989341916e-13*phm4_7 
	- 1.5308084989341916e-13*phm6_7 - 2500.0*phm5_6 + 2500.0*phm5_8)^2 + 
	(-200.0*psm5_7 + 100.0*psm4_7 + 100.0*psm6_7 - 200.0*psm5_7 + 100.0*psm5_6 
	+ 100.0*psm5_8 + 0.25*omm5_7)^2 + (-25.0*(psm5_6-psm5_8) * 
	(phm4_7-phm6_7) + 25.0*(psm4_7-psm6_7) * (phm5_6-phm5_8) - 200.0*phm5_7 
	+ 100.0*phm4_7 + 100.0*phm6_7 - 200.0*phm5_7 + 100.0*phm5_6 + 
	100.0*phm5_8)^2 + (-200.0*omm4_7 + 100.0*omm3_7 + 100.0*omm5_7 - 
	200.0*omm4_7 + 100.0*omm4_6 + 100.0*omm4_8 + 1.5308084989341916e-13*phm3_7 
	- 1.5308084989341916e-13*phm5_7 - 2500.0*phm4_6 + 2500.0*phm4_8)^2 + 
	(-200.0*psm4_7 + 100.0*psm3_7 + 100.0*psm5_7 - 200.0*psm4_7 + 100.0*psm4_6 
	+ 100.0*psm4_8 + 0.25*omm4_7)^2 + (-25.0*(psm4_6-psm4_8) * 
	(phm3_7-phm5_7) + 25.0*(psm3_7-psm5_7) * (phm4_6-phm4_8) - 200.0*phm4_7 
	+ 100.0*phm3_7 + 100.0*phm5_7 - 200.0*phm4_7 + 100.0*phm4_6 + 
	100.0*phm4_8)^2 + (-200.0*omm3_7 + 100.0*omm2_7 + 100.0*omm4_7 - 
	200.0*omm3_7 + 100.0*omm3_6 + 100.0*omm3_8 + 1.5308084989341916e-13*phm2_7 
	- 1.5308084989341916e-13*phm4_7 - 2500.0*phm3_6 + 2500.0*phm3_8)^2 + 
	(-200.0*psm3_7 + 100.0*psm2_7 + 100.0*psm4_7 - 200.0*psm3_7 + 100.0*psm3_6 
	+ 100.0*psm3_8 + 0.25*omm3_7)^2 + (-25.0*(psm3_6-psm3_8) * 
	(phm2_7-phm4_7) + 25.0*(psm2_7-psm4_7) * (phm3_6-phm3_8) - 200.0*phm3_7 
	+ 100.0*phm2_7 + 100.0*phm4_7 - 200.0*phm3_7 + 100.0*phm3_6 + 
	100.0*phm3_8)^2 + (-200.0*omm2_7 + 100.0*omm1_7 + 100.0*omm3_7 - 
	200.0*omm2_7 + 100.0*omm2_6 + 100.0*omm2_8 + 1.5308084989341916e-13*phm1_7 
	- 1.5308084989341916e-13*phm3_7 - 2500.0*phm2_6 + 2500.0*phm2_8)^2 + 
	(-200.0*psm2_7 + 100.0*psm1_7 + 100.0*psm3_7 - 200.0*psm2_7 + 100.0*psm2_6 
	+ 100.0*psm2_8 + 0.25*omm2_7)^2 + (-25.0*(psm2_6-psm2_8) * 
	(phm1_7-phm3_7) + 25.0*(psm1_7-psm3_7) * (phm2_6-phm2_8) - 200.0*phm2_7 
	+ 100.0*phm1_7 + 100.0*phm3_7 - 200.0*phm2_7 + 100.0*phm2_6 + 
	100.0*phm2_8)^2 + (-200.0*omm1_7 + 100.0*om0_7 + 100.0*omm2_7 - 
	200.0*omm1_7 + 100.0*omm1_6 + 100.0*omm1_8 + 1.5308084989341916e-13*ph0_7 - 
	1.5308084989341916e-13*phm2_7 - 2500.0*phm1_6 + 2500.0*phm1_8)^2 + 
	(-200.0*psm1_7 + 100.0*ps0_7 + 100.0*psm2_7 - 200.0*psm1_7 + 100.0*psm1_6 
	+ 100.0*psm1_8 + 0.25*omm1_7)^2 + (-25.0*(psm1_6-psm1_8) * (ph0_7-phm2_7) 
	+ 25.0*(ps0_7-psm2_7) * (phm1_6-phm1_8) - 200.0*phm1_7 + 100.0*ph0_7 + 
	100.0*phm2_7 - 200.0*phm1_7 + 100.0*phm1_6 + 100.0*phm1_8)^2 + 
	(-200.0*om0_7 + 100.0*om1_7 + 100.0*omm1_7 - 200.0*om0_7 + 100.0*om0_6 + 
	100.0*om0_8 + 1.5308084989341916e-13*ph1_7 - 1.5308084989341916e-13*phm1_7 - 
	2500.0*ph0_6 + 2500.0*ph0_8)^2 + (-200.0*ps0_7 + 100.0*ps1_7 + 
	100.0*psm1_7 - 200.0*ps0_7 + 100.0*ps0_6 + 100.0*ps0_8 + 0.25*om0_7)^2 + 
	(-25.0*(ps0_6-ps0_8) * (ph1_7-phm1_7) + 25.0*(ps1_7-psm1_7) * 
	(ph0_6-ph0_8) - 200.0*ph0_7 + 100.0*ph1_7 + 100.0*phm1_7 - 200.0*ph0_7 + 
	100.0*ph0_6 + 100.0*ph0_8)^2 + (-200.0*om1_7 + 100.0*om2_7 + 100.0*om0_7 - 
	200.0*om1_7 + 100.0*om1_6 + 100.0*om1_8 + 1.5308084989341916e-13*ph2_7 - 
	1.5308084989341916e-13*ph0_7 - 2500.0*ph1_6 + 2500.0*ph1_8)^2 + 
	(-200.0*ps1_7 + 100.0*ps2_7 + 100.0*ps0_7 - 200.0*ps1_7 + 100.0*ps1_6 + 
	100.0*ps1_8 + 0.25*om1_7)^2 + (-25.0*(ps1_6-ps1_8) * (ph2_7-ph0_7) + 
	25.0*(ps2_7-ps0_7) * (ph1_6-ph1_8) - 200.0*ph1_7 + 100.0*ph2_7 + 
	100.0*ph0_7 - 200.0*ph1_7 + 100.0*ph1_6 + 100.0*ph1_8)^2 + (-200.0*om2_7 + 
	100.0*om3_7 + 100.0*om1_7 - 200.0*om2_7 + 100.0*om2_6 + 100.0*om2_8 + 
	1.5308084989341916e-13*ph3_7 - 1.5308084989341916e-13*ph1_7 - 2500.0*ph2_6 + 
	2500.0*ph2_8)^2 + (-200.0*ps2_7 + 100.0*ps3_7 + 100.0*ps1_7 - 200.0*ps2_7 
	+ 100.0*ps2_6 + 100.0*ps2_8 + 0.25*om2_7)^2 + (-25.0*(ps2_6-ps2_8) * 
	(ph3_7-ph1_7) + 25.0*(ps3_7-ps1_7) * (ph2_6-ph2_8) - 200.0*ph2_7 + 
	100.0*ph3_7 + 100.0*ph1_7 - 200.0*ph2_7 + 100.0*ph2_6 + 100.0*ph2_8)^2 + 
	(-200.0*om3_7 + 100.0*om4_7 + 100.0*om2_7 - 200.0*om3_7 + 100.0*om3_6 + 
	100.0*om3_8 + 1.5308084989341916e-13*ph4_7 - 1.5308084989341916e-13*ph2_7 - 
	2500.0*ph3_6 + 2500.0*ph3_8)^2 + (-200.0*ps3_7 + 100.0*ps4_7 + 100.0*ps2_7 
	- 200.0*ps3_7 + 100.0*ps3_6 + 100.0*ps3_8 + 0.25*om3_7)^2 + 
	(-25.0*(ps3_6-ps3_8) * (ph4_7-ph2_7) + 25.0*(ps4_7-ps2_7) * 
	(ph3_6-ph3_8) - 200.0*ph3_7 + 100.0*ph4_7 + 100.0*ph2_7 - 200.0*ph3_7 + 
	100.0*ph3_6 + 100.0*ph3_8)^2 + (-200.0*om4_7 + 100.0*om5_7 + 100.0*om3_7 - 
	200.0*om4_7 + 100.0*om4_6 + 100.0*om4_8 + 1.5308084989341916e-13*ph5_7 - 
	1.5308084989341916e-13*ph3_7 - 2500.0*ph4_6 + 2500.0*ph4_8)^2 + 
	(-200.0*ps4_7 + 100.0*ps5_7 + 100.0*ps3_7 - 200.0*ps4_7 + 100.0*ps4_6 + 
	100.0*ps4_8 + 0.25*om4_7)^2 + (-25.0*(ps4_6-ps4_8) * (ph5_7-ph3_7) + 
	25.0*(ps5_7-ps3_7) * (ph4_6-ph4_8) - 200.0*ph4_7 + 100.0*ph5_7 + 
	100.0*ph3_7 - 200.0*ph4_7 + 100.0*ph4_6 + 100.0*ph4_8)^2 + (-200.0*om5_7 + 
	100.0*om6_7 + 100.0*om4_7 - 200.0*om5_7 + 100.0*om5_6 + 100.0*om5_8 + 
	1.5308084989341916e-13*ph6_7 - 1.5308084989341916e-13*ph4_7 - 2500.0*ph5_6 + 
	2500.0*ph5_8)^2 + (-200.0*ps5_7 + 100.0*ps6_7 + 100.0*ps4_7 - 200.0*ps5_7 
	+ 100.0*ps5_6 + 100.0*ps5_8 + 0.25*om5_7)^2 + (-25.0*(ps5_6-ps5_8) * 
	(ph6_7-ph4_7) + 25.0*(ps6_7-ps4_7) * (ph5_6-ph5_8) - 200.0*ph5_7 + 
	100.0*ph6_7 + 100.0*ph4_7 - 200.0*ph5_7 + 100.0*ph5_6 + 100.0*ph5_8)^2 + 
	(-200.0*om6_7 + 100.0*om7_7 + 100.0*om5_7 - 200.0*om6_7 + 100.0*om6_6 + 
	100.0*om6_8 + 1.5308084989341916e-13*ph7_7 - 1.5308084989341916e-13*ph5_7 - 
	2500.0*ph6_6 + 2500.0*ph6_8)^2 + (-200.0*ps6_7 + 100.0*ps7_7 + 100.0*ps5_7 
	- 200.0*ps6_7 + 100.0*ps6_6 + 100.0*ps6_8 + 0.25*om6_7)^2 + 
	(-25.0*(ps6_6-ps6_8) * (ph7_7-ph5_7) + 25.0*(ps7_7-ps5_7) * 
	(ph6_6-ph6_8) - 200.0*ph6_7 + 100.0*ph7_7 + 100.0*ph5_7 - 200.0*ph6_7 + 
	100.0*ph6_6 + 100.0*ph6_8)^2 + (-200.0*om7_7 + 100.0*om8_7 + 100.0*om6_7 - 
	200.0*om7_7 + 100.0*om7_6 + 100.0*om7_8 + 1.5308084989341916e-13*ph8_7 - 
	1.5308084989341916e-13*ph6_7 - 2500.0*ph7_6 + 2500.0*ph7_8)^2 + 
	(-200.0*ps7_7 + 100.0*ps8_7 + 100.0*ps6_7 - 200.0*ps7_7 + 100.0*ps7_6 + 
	100.0*ps7_8 + 0.25*om7_7)^2 + (-25.0*(ps7_6-ps7_8) * (ph8_7-ph6_7) + 
	25.0*(ps8_7-ps6_7) * (ph7_6-ph7_8) - 200.0*ph7_7 + 100.0*ph8_7 + 
	100.0*ph6_7 - 200.0*ph7_7 + 100.0*ph7_6 + 100.0*ph7_8)^2 + (-200.0*om8_7 + 
	100.0*om9_7 + 100.0*om7_7 - 200.0*om8_7 + 100.0*om8_6 + 100.0*om8_8 + 
	1.5308084989341916e-13*ph9_7 - 1.5308084989341916e-13*ph7_7 - 2500.0*ph8_6 + 
	2500.0*ph8_8)^2 + (-200.0*ps8_7 + 100.0*ps9_7 + 100.0*ps7_7 - 200.0*ps8_7 
	+ 100.0*ps8_6 + 100.0*ps8_8 + 0.25*om8_7)^2 + (-25.0*(ps8_6-ps8_8) * 
	(ph9_7-ph7_7) + 25.0*(ps9_7-ps7_7) * (ph8_6-ph8_8) - 200.0*ph8_7 + 
	100.0*ph9_7 + 100.0*ph7_7 - 200.0*ph8_7 + 100.0*ph8_6 + 100.0*ph8_8)^2 + 
	(-200.0*om9_7 + 100.0*om10_7 + 100.0*om8_7 - 200.0*om9_7 + 100.0*om9_6 + 
	100.0*om9_8 + 1.5308084989341916e-13*ph10_7 - 1.5308084989341916e-13*ph8_7 - 
	2500.0*ph9_6 + 2500.0*ph9_8)^2 + (-200.0*ps9_7 + 100.0*ps10_7 + 
	100.0*ps8_7 - 200.0*ps9_7 + 100.0*ps9_6 + 100.0*ps9_8 + 0.25*om9_7)^2 + 
	(-25.0*(ps9_6-ps9_8) * (ph10_7-ph8_7) + 25.0*(ps10_7-ps8_7) * 
	(ph9_6-ph9_8) - 200.0*ph9_7 + 100.0*ph10_7 + 100.0*ph8_7 - 200.0*ph9_7 + 
	100.0*ph9_6 + 100.0*ph9_8)^2 + (-200.0*omm9_6 + 100.0*omm8_6 + 
	100.0*omm10_6 - 200.0*omm9_6 + 100.0*omm9_5 + 100.0*omm9_7 + 
	1.5308084989341916e-13*phm8_6 - 1.5308084989341916e-13*phm10_6 - 
	2500.0*phm9_5 + 2500.0*phm9_7)^2 + (-200.0*psm9_6 + 100.0*psm8_6 + 
	100.0*psm10_6 - 200.0*psm9_6 + 100.0*psm9_5 + 100.0*psm9_7 + 
	0.25*omm9_6)^2 + (-25.0*(psm9_5-psm9_7) * (phm8_6-phm10_6) + 
	25.0*(psm8_6-psm10_6) * (phm9_5-phm9_7) - 200.0*phm9_6 + 100.0*phm8_6 + 
	100.0*phm10_6 - 200.0*phm9_6 + 100.0*phm9_5 + 100.0*phm9_7)^2 + 
	(-200.0*omm8_6 + 100.0*omm7_6 + 100.0*omm9_6 - 200.0*omm8_6 + 100.0*omm8_5 
	+ 100.0*omm8_7 + 1.5308084989341916e-13*phm7_6 - 
	1.5308084989341916e-13*phm9_6 - 2500.0*phm8_5 + 2500.0*phm8_7)^2 + 
	(-200.0*psm8_6 + 100.0*psm7_6 + 100.0*psm9_6 - 200.0*psm8_6 + 100.0*psm8_5 
	+ 100.0*psm8_7 + 0.25*omm8_6)^2 + (-25.0*(psm8_5-psm8_7) * 
	(phm7_6-phm9_6) + 25.0*(psm7_6-psm9_6) * (phm8_5-phm8_7) - 200.0*phm8_6 
	+ 100.0*phm7_6 + 100.0*phm9_6 - 200.0*phm8_6 + 100.0*phm8_5 + 
	100.0*phm8_7)^2 + (-200.0*omm7_6 + 100.0*omm6_6 + 100.0*omm8_6 - 
	200.0*omm7_6 + 100.0*omm7_5 + 100.0*omm7_7 + 1.5308084989341916e-13*phm6_6 
	- 1.5308084989341916e-13*phm8_6 - 2500.0*phm7_5 + 2500.0*phm7_7)^2 + 
	(-200.0*psm7_6 + 100.0*psm6_6 + 100.0*psm8_6 - 200.0*psm7_6 + 100.0*psm7_5 
	+ 100.0*psm7_7 + 0.25*omm7_6)^2 + (-25.0*(psm7_5-psm7_7) * 
	(phm6_6-phm8_6) + 25.0*(psm6_6-psm8_6) * (phm7_5-phm7_7) - 200.0*phm7_6 
	+ 100.0*phm6_6 + 100.0*phm8_6 - 200.0*phm7_6 + 100.0*phm7_5 + 
	100.0*phm7_7)^2 + (-200.0*omm6_6 + 100.0*omm5_6 + 100.0*omm7_6 - 
	200.0*omm6_6 + 100.0*omm6_5 + 100.0*omm6_7 + 1.5308084989341916e-13*phm5_6 
	- 1.5308084989341916e-13*phm7_6 - 2500.0*phm6_5 + 2500.0*phm6_7)^2 + 
	(-200.0*psm6_6 + 100.0*psm5_6 + 100.0*psm7_6 - 200.0*psm6_6 + 100.0*psm6_5 
	+ 100.0*psm6_7 + 0.25*omm6_6)^2 + (-25.0*(psm6_5-psm6_7) * 
	(phm5_6-phm7_6) + 25.0*(psm5_6-psm7_6) * (phm6_5-phm6_7) - 200.0*phm6_6 
	+ 100.0*phm5_6 + 100.0*phm7_6 - 200.0*phm6_6 + 100.0*phm6_5 + 
	100.0*phm6_7)^2 + (-200.0*omm5_6 + 100.0*omm4_6 + 100.0*omm6_6 - 
	200.0*omm5_6 + 100.0*omm5_5 + 100.0*omm5_7 + 1.5308084989341916e-13*phm4_6 
	- 1.5308084989341916e-13*phm6_6 - 2500.0*phm5_5 + 2500.0*phm5_7)^2 + 
	(-200.0*psm5_6 + 100.0*psm4_6 + 100.0*psm6_6 - 200.0*psm5_6 + 100.0*psm5_5 
	+ 100.0*psm5_7 + 0.25*omm5_6)^2 + (-25.0*(psm5_5-psm5_7) * 
	(phm4_6-phm6_6) + 25.0*(psm4_6-psm6_6) * (phm5_5-phm5_7) - 200.0*phm5_6 
	+ 100.0*phm4_6 + 100.0*phm6_6 - 200.0*phm5_6 + 100.0*phm5_5 + 
	100.0*phm5_7)^2 + (-200.0*omm4_6 + 100.0*omm3_6 + 100.0*omm5_6 - 
	200.0*omm4_6 + 100.0*omm4_5 + 100.0*omm4_7 + 1.5308084989341916e-13*phm3_6 
	- 1.5308084989341916e-13*phm5_6 - 2500.0*phm4_5 + 2500.0*phm4_7)^2 + 
	(-200.0*psm4_6 + 100.0*psm3_6 + 100.0*psm5_6 - 200.0*psm4_6 + 100.0*psm4_5 
	+ 100.0*psm4_7 + 0.25*omm4_6)^2 + (-25.0*(psm4_5-psm4_7) * 
	(phm3_6-phm5_6) + 25.0*(psm3_6-psm5_6) * (phm4_5-phm4_7) - 200.0*phm4_6 
	+ 100.0*phm3_6 + 100.0*phm5_6 - 200.0*phm4_6 + 100.0*phm4_5 + 
	100.0*phm4_7)^2 + (-200.0*omm3_6 + 100.0*omm2_6 + 100.0*omm4_6 - 
	200.0*omm3_6 + 100.0*omm3_5 + 100.0*omm3_7 + 1.5308084989341916e-13*phm2_6 
	- 1.5308084989341916e-13*phm4_6 - 2500.0*phm3_5 + 2500.0*phm3_7)^2 + 
	(-200.0*psm3_6 + 100.0*psm2_6 + 100.0*psm4_6 - 200.0*psm3_6 + 100.0*psm3_5 
	+ 100.0*psm3_7 + 0.25*omm3_6)^2 + (-25.0*(psm3_5-psm3_7) * 
	(phm2_6-phm4_6) + 25.0*(psm2_6-psm4_6) * (phm3_5-phm3_7) - 200.0*phm3_6 
	+ 100.0*phm2_6 + 100.0*phm4_6 - 200.0*phm3_6 + 100.0*phm3_5 + 
	100.0*phm3_7)^2 + (-200.0*omm2_6 + 100.0*omm1_6 + 100.0*omm3_6 - 
	200.0*omm2_6 + 100.0*omm2_5 + 100.0*omm2_7 + 1.5308084989341916e-13*phm1_6 
	- 1.5308084989341916e-13*phm3_6 - 2500.0*phm2_5 + 2500.0*phm2_7)^2 + 
	(-200.0*psm2_6 + 100.0*psm1_6 + 100.0*psm3_6 - 200.0*psm2_6 + 100.0*psm2_5 
	+ 100.0*psm2_7 + 0.25*omm2_6)^2 + (-25.0*(psm2_5-psm2_7) * 
	(phm1_6-phm3_6) + 25.0*(psm1_6-psm3_6) * (phm2_5-phm2_7) - 200.0*phm2_6 
	+ 100.0*phm1_6 + 100.0*phm3_6 - 200.0*phm2_6 + 100.0*phm2_5 + 
	100.0*phm2_7)^2 + (-200.0*omm1_6 + 100.0*om0_6 + 100.0*omm2_6 - 
	200.0*omm1_6 + 100.0*omm1_5 + 100.0*omm1_7 + 1.5308084989341916e-13*ph0_6 - 
	1.5308084989341916e-13*phm2_6 - 2500.0*phm1_5 + 2500.0*phm1_7)^2 + 
	(-200.0*psm1_6 + 100.0*ps0_6 + 100.0*psm2_6 - 200.0*psm1_6 + 100.0*psm1_5 
	+ 100.0*psm1_7 + 0.25*omm1_6)^2 + (-25.0*(psm1_5-psm1_7) * (ph0_6-phm2_6) 
	+ 25.0*(ps0_6-psm2_6) * (phm1_5-phm1_7) - 200.0*phm1_6 + 100.0*ph0_6 + 
	100.0*phm2_6 - 200.0*phm1_6 + 100.0*phm1_5 + 100.0*phm1_7)^2 + 
	(-200.0*om0_6 + 100.0*om1_6 + 100.0*omm1_6 - 200.0*om0_6 + 100.0*om0_5 + 
	100.0*om0_7 + 1.5308084989341916e-13*ph1_6 - 1.5308084989341916e-13*phm1_6 - 
	2500.0*ph0_5 + 2500.0*ph0_7)^2 + (-200.0*ps0_6 + 100.0*ps1_6 + 
	100.0*psm1_6 - 200.0*ps0_6 + 100.0*ps0_5 + 100.0*ps0_7 + 0.25*om0_6)^2 + 
	(-25.0*(ps0_5-ps0_7) * (ph1_6-phm1_6) + 25.0*(ps1_6-psm1_6) * 
	(ph0_5-ph0_7) - 200.0*ph0_6 + 100.0*ph1_6 + 100.0*phm1_6 - 200.0*ph0_6 + 
	100.0*ph0_5 + 100.0*ph0_7)^2 + (-200.0*om1_6 + 100.0*om2_6 + 100.0*om0_6 - 
	200.0*om1_6 + 100.0*om1_5 + 100.0*om1_7 + 1.5308084989341916e-13*ph2_6 - 
	1.5308084989341916e-13*ph0_6 - 2500.0*ph1_5 + 2500.0*ph1_7)^2 + 
	(-200.0*ps1_6 + 100.0*ps2_6 + 100.0*ps0_6 - 200.0*ps1_6 + 100.0*ps1_5 + 
	100.0*ps1_7 + 0.25*om1_6)^2 + (-25.0*(ps1_5-ps1_7) * (ph2_6-ph0_6) + 
	25.0*(ps2_6-ps0_6) * (ph1_5-ph1_7) - 200.0*ph1_6 + 100.0*ph2_6 + 
	100.0*ph0_6 - 200.0*ph1_6 + 100.0*ph1_5 + 100.0*ph1_7)^2 + (-200.0*om2_6 + 
	100.0*om3_6 + 100.0*om1_6 - 200.0*om2_6 + 100.0*om2_5 + 100.0*om2_7 + 
	1.5308084989341916e-13*ph3_6 - 1.5308084989341916e-13*ph1_6 - 2500.0*ph2_5 + 
	2500.0*ph2_7)^2 + (-200.0*ps2_6 + 100.0*ps3_6 + 100.0*ps1_6 - 200.0*ps2_6 
	+ 100.0*ps2_5 + 100.0*ps2_7 + 0.25*om2_6)^2 + (-25.0*(ps2_5-ps2_7) * 
	(ph3_6-ph1_6) + 25.0*(ps3_6-ps1_6) * (ph2_5-ph2_7) - 200.0*ph2_6 + 
	100.0*ph3_6 + 100.0*ph1_6 - 200.0*ph2_6 + 100.0*ph2_5 + 100.0*ph2_7)^2 + 
	(-200.0*om3_6 + 100.0*om4_6 + 100.0*om2_6 - 200.0*om3_6 + 100.0*om3_5 + 
	100.0*om3_7 + 1.5308084989341916e-13*ph4_6 - 1.5308084989341916e-13*ph2_6 - 
	2500.0*ph3_5 + 2500.0*ph3_7)^2 + (-200.0*ps3_6 + 100.0*ps4_6 + 100.0*ps2_6 
	- 200.0*ps3_6 + 100.0*ps3_5 + 100.0*ps3_7 + 0.25*om3_6)^2 + 
	(-25.0*(ps3_5-ps3_7) * (ph4_6-ph2_6) + 25.0*(ps4_6-ps2_6) * 
	(ph3_5-ph3_7) - 200.0*ph3_6 + 100.0*ph4_6 + 100.0*ph2_6 - 200.0*ph3_6 + 
	100.0*ph3_5 + 100.0*ph3_7)^2 + (-200.0*om4_6 + 100.0*om5_6 + 100.0*om3_6 - 
	200.0*om4_6 + 100.0*om4_5 + 100.0*om4_7 + 1.5308084989341916e-13*ph5_6 - 
	1.5308084989341916e-13*ph3_6 - 2500.0*ph4_5 + 2500.0*ph4_7)^2 + 
	(-200.0*ps4_6 + 100.0*ps5_6 + 100.0*ps3_6 - 200.0*ps4_6 + 100.0*ps4_5 + 
	100.0*ps4_7 + 0.25*om4_6)^2 + (-25.0*(ps4_5-ps4_7) * (ph5_6-ph3_6) + 
	25.0*(ps5_6-ps3_6) * (ph4_5-ph4_7) - 200.0*ph4_6 + 100.0*ph5_6 + 
	100.0*ph3_6 - 200.0*ph4_6 + 100.0*ph4_5 + 100.0*ph4_7)^2 + (-200.0*om5_6 + 
	100.0*om6_6 + 100.0*om4_6 - 200.0*om5_6 + 100.0*om5_5 + 100.0*om5_7 + 
	1.5308084989341916e-13*ph6_6 - 1.5308084989341916e-13*ph4_6 - 2500.0*ph5_5 + 
	2500.0*ph5_7)^2 + (-200.0*ps5_6 + 100.0*ps6_6 + 100.0*ps4_6 - 200.0*ps5_6 
	+ 100.0*ps5_5 + 100.0*ps5_7 + 0.25*om5_6)^2 + (-25.0*(ps5_5-ps5_7) * 
	(ph6_6-ph4_6) + 25.0*(ps6_6-ps4_6) * (ph5_5-ph5_7) - 200.0*ph5_6 + 
	100.0*ph6_6 + 100.0*ph4_6 - 200.0*ph5_6 + 100.0*ph5_5 + 100.0*ph5_7)^2 + 
	(-200.0*om6_6 + 100.0*om7_6 + 100.0*om5_6 - 200.0*om6_6 + 100.0*om6_5 + 
	100.0*om6_7 + 1.5308084989341916e-13*ph7_6 - 1.5308084989341916e-13*ph5_6 - 
	2500.0*ph6_5 + 2500.0*ph6_7)^2 + (-200.0*ps6_6 + 100.0*ps7_6 + 100.0*ps5_6 
	- 200.0*ps6_6 + 100.0*ps6_5 + 100.0*ps6_7 + 0.25*om6_6)^2 + 
	(-25.0*(ps6_5-ps6_7) * (ph7_6-ph5_6) + 25.0*(ps7_6-ps5_6) * 
	(ph6_5-ph6_7) - 200.0*ph6_6 + 100.0*ph7_6 + 100.0*ph5_6 - 200.0*ph6_6 + 
	100.0*ph6_5 + 100.0*ph6_7)^2 + (-200.0*om7_6 + 100.0*om8_6 + 100.0*om6_6 - 
	200.0*om7_6 + 100.0*om7_5 + 100.0*om7_7 + 1.5308084989341916e-13*ph8_6 - 
	1.5308084989341916e-13*ph6_6 - 2500.0*ph7_5 + 2500.0*ph7_7)^2 + 
	(-200.0*ps7_6 + 100.0*ps8_6 + 100.0*ps6_6 - 200.0*ps7_6 + 100.0*ps7_5 + 
	100.0*ps7_7 + 0.25*om7_6)^2 + (-25.0*(ps7_5-ps7_7) * (ph8_6-ph6_6) + 
	25.0*(ps8_6-ps6_6) * (ph7_5-ph7_7) - 200.0*ph7_6 + 100.0*ph8_6 + 
	100.0*ph6_6 - 200.0*ph7_6 + 100.0*ph7_5 + 100.0*ph7_7)^2 + (-200.0*om8_6 + 
	100.0*om9_6 + 100.0*om7_6 - 200.0*om8_6 + 100.0*om8_5 + 100.0*om8_7 + 
	1.5308084989341916e-13*ph9_6 - 1.5308084989341916e-13*ph7_6 - 2500.0*ph8_5 + 
	2500.0*ph8_7)^2 + (-200.0*ps8_6 + 100.0*ps9_6 + 100.0*ps7_6 - 200.0*ps8_6 
	+ 100.0*ps8_5 + 100.0*ps8_7 + 0.25*om8_6)^2 + (-25.0*(ps8_5-ps8_7) * 
	(ph9_6-ph7_6) + 25.0*(ps9_6-ps7_6) * (ph8_5-ph8_7) - 200.0*ph8_6 + 
	100.0*ph9_6 + 100.0*ph7_6 - 200.0*ph8_6 + 100.0*ph8_5 + 100.0*ph8_7)^2 + 
	(-200.0*om9_6 + 100.0*om10_6 + 100.0*om8_6 - 200.0*om9_6 + 100.0*om9_5 + 
	100.0*om9_7 + 1.5308084989341916e-13*ph10_6 - 1.5308084989341916e-13*ph8_6 - 
	2500.0*ph9_5 + 2500.0*ph9_7)^2 + (-200.0*ps9_6 + 100.0*ps10_6 + 
	100.0*ps8_6 - 200.0*ps9_6 + 100.0*ps9_5 + 100.0*ps9_7 + 0.25*om9_6)^2 + 
	(-25.0*(ps9_5-ps9_7) * (ph10_6-ph8_6) + 25.0*(ps10_6-ps8_6) * 
	(ph9_5-ph9_7) - 200.0*ph9_6 + 100.0*ph10_6 + 100.0*ph8_6 - 200.0*ph9_6 + 
	100.0*ph9_5 + 100.0*ph9_7)^2 + (-200.0*omm9_5 + 100.0*omm8_5 + 
	100.0*omm10_5 - 200.0*omm9_5 + 100.0*omm9_4 + 100.0*omm9_6 + 
	1.5308084989341916e-13*phm8_5 - 1.5308084989341916e-13*phm10_5 - 
	2500.0*phm9_4 + 2500.0*phm9_6)^2 + (-200.0*psm9_5 + 100.0*psm8_5 + 
	100.0*psm10_5 - 200.0*psm9_5 + 100.0*psm9_4 + 100.0*psm9_6 + 
	0.25*omm9_5)^2 + (-25.0*(psm9_4-psm9_6) * (phm8_5-phm10_5) + 
	25.0*(psm8_5-psm10_5) * (phm9_4-phm9_6) - 200.0*phm9_5 + 100.0*phm8_5 + 
	100.0*phm10_5 - 200.0*phm9_5 + 100.0*phm9_4 + 100.0*phm9_6)^2 + 
	(-200.0*omm8_5 + 100.0*omm7_5 + 100.0*omm9_5 - 200.0*omm8_5 + 100.0*omm8_4 
	+ 100.0*omm8_6 + 1.5308084989341916e-13*phm7_5 - 
	1.5308084989341916e-13*phm9_5 - 2500.0*phm8_4 + 2500.0*phm8_6)^2 + 
	(-200.0*psm8_5 + 100.0*psm7_5 + 100.0*psm9_5 - 200.0*psm8_5 + 100.0*psm8_4 
	+ 100.0*psm8_6 + 0.25*omm8_5)^2 + (-25.0*(psm8_4-psm8_6) * 
	(phm7_5-phm9_5) + 25.0*(psm7_5-psm9_5) * (phm8_4-phm8_6) - 200.0*phm8_5 
	+ 100.0*phm7_5 + 100.0*phm9_5 - 200.0*phm8_5 + 100.0*phm8_4 + 
	100.0*phm8_6)^2 + (-200.0*omm7_5 + 100.0*omm6_5 + 100.0*omm8_5 - 
	200.0*omm7_5 + 100.0*omm7_4 + 100.0*omm7_6 + 1.5308084989341916e-13*phm6_5 
	- 1.5308084989341916e-13*phm8_5 - 2500.0*phm7_4 + 2500.0*phm7_6)^2 + 
	(-200.0*psm7_5 + 100.0*psm6_5 + 100.0*psm8_5 - 200.0*psm7_5 + 100.0*psm7_4 
	+ 100.0*psm7_6 + 0.25*omm7_5)^2 + (-25.0*(psm7_4-psm7_6) * 
	(phm6_5-phm8_5) + 25.0*(psm6_5-psm8_5) * (phm7_4-phm7_6) - 200.0*phm7_5 
	+ 100.0*phm6_5 + 100.0*phm8_5 - 200.0*phm7_5 + 100.0*phm7_4 + 
	100.0*phm7_6)^2 + (-200.0*omm6_5 + 100.0*omm5_5 + 100.0*omm7_5 - 
	200.0*omm6_5 + 100.0*omm6_4 + 100.0*omm6_6 + 1.5308084989341916e-13*phm5_5 
	- 1.5308084989341916e-13*phm7_5 - 2500.0*phm6_4 + 2500.0*phm6_6)^2 + 
	(-200.0*psm6_5 + 100.0*psm5_5 + 100.0*psm7_5 - 200.0*psm6_5 + 100.0*psm6_4 
	+ 100.0*psm6_6 + 0.25*omm6_5)^2 + (-25.0*(psm6_4-psm6_6) * 
	(phm5_5-phm7_5) + 25.0*(psm5_5-psm7_5) * (phm6_4-phm6_6) - 200.0*phm6_5 
	+ 100.0*phm5_5 + 100.0*phm7_5 - 200.0*phm6_5 + 100.0*phm6_4 + 
	100.0*phm6_6)^2 + (-200.0*omm5_5 + 100.0*omm4_5 + 100.0*omm6_5 - 
	200.0*omm5_5 + 100.0*omm5_4 + 100.0*omm5_6 + 1.5308084989341916e-13*phm4_5 
	- 1.5308084989341916e-13*phm6_5 - 2500.0*phm5_4 + 2500.0*phm5_6)^2 + 
	(-200.0*psm5_5 + 100.0*psm4_5 + 100.0*psm6_5 - 200.0*psm5_5 + 100.0*psm5_4 
	+ 100.0*psm5_6 + 0.25*omm5_5)^2 + (-25.0*(psm5_4-psm5_6) * 
	(phm4_5-phm6_5) + 25.0*(psm4_5-psm6_5) * (phm5_4-phm5_6) - 200.0*phm5_5 
	+ 100.0*phm4_5 + 100.0*phm6_5 - 200.0*phm5_5 + 100.0*phm5_4 + 
	100.0*phm5_6)^2 + (-200.0*omm4_5 + 100.0*omm3_5 + 100.0*omm5_5 - 
	200.0*omm4_5 + 100.0*omm4_4 + 100.0*omm4_6 + 1.5308084989341916e-13*phm3_5 
	- 1.5308084989341916e-13*phm5_5 - 2500.0*phm4_4 + 2500.0*phm4_6)^2 + 
	(-200.0*psm4_5 + 100.0*psm3_5 + 100.0*psm5_5 - 200.0*psm4_5 + 100.0*psm4_4 
	+ 100.0*psm4_6 + 0.25*omm4_5)^2 + (-25.0*(psm4_4-psm4_6) * 
	(phm3_5-phm5_5) + 25.0*(psm3_5-psm5_5) * (phm4_4-phm4_6) - 200.0*phm4_5 
	+ 100.0*phm3_5 + 100.0*phm5_5 - 200.0*phm4_5 + 100.0*phm4_4 + 
	100.0*phm4_6)^2 + (-200.0*omm3_5 + 100.0*omm2_5 + 100.0*omm4_5 - 
	200.0*omm3_5 + 100.0*omm3_4 + 100.0*omm3_6 + 1.5308084989341916e-13*phm2_5 
	- 1.5308084989341916e-13*phm4_5 - 2500.0*phm3_4 + 2500.0*phm3_6)^2 + 
	(-200.0*psm3_5 + 100.0*psm2_5 + 100.0*psm4_5 - 200.0*psm3_5 + 100.0*psm3_4 
	+ 100.0*psm3_6 + 0.25*omm3_5)^2 + (-25.0*(psm3_4-psm3_6) * 
	(phm2_5-phm4_5) + 25.0*(psm2_5-psm4_5) * (phm3_4-phm3_6) - 200.0*phm3_5 
	+ 100.0*phm2_5 + 100.0*phm4_5 - 200.0*phm3_5 + 100.0*phm3_4 + 
	100.0*phm3_6)^2 + (-200.0*omm2_5 + 100.0*omm1_5 + 100.0*omm3_5 - 
	200.0*omm2_5 + 100.0*omm2_4 + 100.0*omm2_6 + 1.5308084989341916e-13*phm1_5 
	- 1.5308084989341916e-13*phm3_5 - 2500.0*phm2_4 + 2500.0*phm2_6)^2 + 
	(-200.0*psm2_5 + 100.0*psm1_5 + 100.0*psm3_5 - 200.0*psm2_5 + 100.0*psm2_4 
	+ 100.0*psm2_6 + 0.25*omm2_5)^2 + (-25.0*(psm2_4-psm2_6) * 
	(phm1_5-phm3_5) + 25.0*(psm1_5-psm3_5) * (phm2_4-phm2_6) - 200.0*phm2_5 
	+ 100.0*phm1_5 + 100.0*phm3_5 - 200.0*phm2_5 + 100.0*phm2_4 + 
	100.0*phm2_6)^2 + (-200.0*omm1_5 + 100.0*om0_5 + 100.0*omm2_5 - 
	200.0*omm1_5 + 100.0*omm1_4 + 100.0*omm1_6 + 1.5308084989341916e-13*ph0_5 - 
	1.5308084989341916e-13*phm2_5 - 2500.0*phm1_4 + 2500.0*phm1_6)^2 + 
	(-200.0*psm1_5 + 100.0*ps0_5 + 100.0*psm2_5 - 200.0*psm1_5 + 100.0*psm1_4 
	+ 100.0*psm1_6 + 0.25*omm1_5)^2 + (-25.0*(psm1_4-psm1_6) * (ph0_5-phm2_5) 
	+ 25.0*(ps0_5-psm2_5) * (phm1_4-phm1_6) - 200.0*phm1_5 + 100.0*ph0_5 + 
	100.0*phm2_5 - 200.0*phm1_5 + 100.0*phm1_4 + 100.0*phm1_6)^2 + 
	(-200.0*om0_5 + 100.0*om1_5 + 100.0*omm1_5 - 200.0*om0_5 + 100.0*om0_4 + 
	100.0*om0_6 + 1.5308084989341916e-13*ph1_5 - 1.5308084989341916e-13*phm1_5 - 
	2500.0*ph0_4 + 2500.0*ph0_6)^2 + (-200.0*ps0_5 + 100.0*ps1_5 + 
	100.0*psm1_5 - 200.0*ps0_5 + 100.0*ps0_4 + 100.0*ps0_6 + 0.25*om0_5)^2 + 
	(-25.0*(ps0_4-ps0_6) * (ph1_5-phm1_5) + 25.0*(ps1_5-psm1_5) * 
	(ph0_4-ph0_6) - 200.0*ph0_5 + 100.0*ph1_5 + 100.0*phm1_5 - 200.0*ph0_5 + 
	100.0*ph0_4 + 100.0*ph0_6)^2 + (-200.0*om1_5 + 100.0*om2_5 + 100.0*om0_5 - 
	200.0*om1_5 + 100.0*om1_4 + 100.0*om1_6 + 1.5308084989341916e-13*ph2_5 - 
	1.5308084989341916e-13*ph0_5 - 2500.0*ph1_4 + 2500.0*ph1_6)^2 + 
	(-200.0*ps1_5 + 100.0*ps2_5 + 100.0*ps0_5 - 200.0*ps1_5 + 100.0*ps1_4 + 
	100.0*ps1_6 + 0.25*om1_5)^2 + (-25.0*(ps1_4-ps1_6) * (ph2_5-ph0_5) + 
	25.0*(ps2_5-ps0_5) * (ph1_4-ph1_6) - 200.0*ph1_5 + 100.0*ph2_5 + 
	100.0*ph0_5 - 200.0*ph1_5 + 100.0*ph1_4 + 100.0*ph1_6)^2 + (-200.0*om2_5 + 
	100.0*om3_5 + 100.0*om1_5 - 200.0*om2_5 + 100.0*om2_4 + 100.0*om2_6 + 
	1.5308084989341916e-13*ph3_5 - 1.5308084989341916e-13*ph1_5 - 2500.0*ph2_4 + 
	2500.0*ph2_6)^2 + (-200.0*ps2_5 + 100.0*ps3_5 + 100.0*ps1_5 - 200.0*ps2_5 
	+ 100.0*ps2_4 + 100.0*ps2_6 + 0.25*om2_5)^2 + (-25.0*(ps2_4-ps2_6) * 
	(ph3_5-ph1_5) + 25.0*(ps3_5-ps1_5) * (ph2_4-ph2_6) - 200.0*ph2_5 + 
	100.0*ph3_5 + 100.0*ph1_5 - 200.0*ph2_5 + 100.0*ph2_4 + 100.0*ph2_6)^2 + 
	(-200.0*om3_5 + 100.0*om4_5 + 100.0*om2_5 - 200.0*om3_5 + 100.0*om3_4 + 
	100.0*om3_6 + 1.5308084989341916e-13*ph4_5 - 1.5308084989341916e-13*ph2_5 - 
	2500.0*ph3_4 + 2500.0*ph3_6)^2 + (-200.0*ps3_5 + 100.0*ps4_5 + 100.0*ps2_5 
	- 200.0*ps3_5 + 100.0*ps3_4 + 100.0*ps3_6 + 0.25*om3_5)^2 + 
	(-25.0*(ps3_4-ps3_6) * (ph4_5-ph2_5) + 25.0*(ps4_5-ps2_5) * 
	(ph3_4-ph3_6) - 200.0*ph3_5 + 100.0*ph4_5 + 100.0*ph2_5 - 200.0*ph3_5 + 
	100.0*ph3_4 + 100.0*ph3_6)^2 + (-200.0*om4_5 + 100.0*om5_5 + 100.0*om3_5 - 
	200.0*om4_5 + 100.0*om4_4 + 100.0*om4_6 + 1.5308084989341916e-13*ph5_5 - 
	1.5308084989341916e-13*ph3_5 - 2500.0*ph4_4 + 2500.0*ph4_6)^2 + 
	(-200.0*ps4_5 + 100.0*ps5_5 + 100.0*ps3_5 - 200.0*ps4_5 + 100.0*ps4_4 + 
	100.0*ps4_6 + 0.25*om4_5)^2 + (-25.0*(ps4_4-ps4_6) * (ph5_5-ph3_5) + 
	25.0*(ps5_5-ps3_5) * (ph4_4-ph4_6) - 200.0*ph4_5 + 100.0*ph5_5 + 
	100.0*ph3_5 - 200.0*ph4_5 + 100.0*ph4_4 + 100.0*ph4_6)^2 + (-200.0*om5_5 + 
	100.0*om6_5 + 100.0*om4_5 - 200.0*om5_5 + 100.0*om5_4 + 100.0*om5_6 + 
	1.5308084989341916e-13*ph6_5 - 1.5308084989341916e-13*ph4_5 - 2500.0*ph5_4 + 
	2500.0*ph5_6)^2 + (-200.0*ps5_5 + 100.0*ps6_5 + 100.0*ps4_5 - 200.0*ps5_5 
	+ 100.0*ps5_4 + 100.0*ps5_6 + 0.25*om5_5)^2 + (-25.0*(ps5_4-ps5_6) * 
	(ph6_5-ph4_5) + 25.0*(ps6_5-ps4_5) * (ph5_4-ph5_6) - 200.0*ph5_5 + 
	100.0*ph6_5 + 100.0*ph4_5 - 200.0*ph5_5 + 100.0*ph5_4 + 100.0*ph5_6)^2 + 
	(-200.0*om6_5 + 100.0*om7_5 + 100.0*om5_5 - 200.0*om6_5 + 100.0*om6_4 + 
	100.0*om6_6 + 1.5308084989341916e-13*ph7_5 - 1.5308084989341916e-13*ph5_5 - 
	2500.0*ph6_4 + 2500.0*ph6_6)^2 + (-200.0*ps6_5 + 100.0*ps7_5 + 100.0*ps5_5 
	- 200.0*ps6_5 + 100.0*ps6_4 + 100.0*ps6_6 + 0.25*om6_5)^2 + 
	(-25.0*(ps6_4-ps6_6) * (ph7_5-ph5_5) + 25.0*(ps7_5-ps5_5) * 
	(ph6_4-ph6_6) - 200.0*ph6_5 + 100.0*ph7_5 + 100.0*ph5_5 - 200.0*ph6_5 + 
	100.0*ph6_4 + 100.0*ph6_6)^2 + (-200.0*om7_5 + 100.0*om8_5 + 100.0*om6_5 - 
	200.0*om7_5 + 100.0*om7_4 + 100.0*om7_6 + 1.5308084989341916e-13*ph8_5 - 
	1.5308084989341916e-13*ph6_5 - 2500.0*ph7_4 + 2500.0*ph7_6)^2 + 
	(-200.0*ps7_5 + 100.0*ps8_5 + 100.0*ps6_5 - 200.0*ps7_5 + 100.0*ps7_4 + 
	100.0*ps7_6 + 0.25*om7_5)^2 + (-25.0*(ps7_4-ps7_6) * (ph8_5-ph6_5) + 
	25.0*(ps8_5-ps6_5) * (ph7_4-ph7_6) - 200.0*ph7_5 + 100.0*ph8_5 + 
	100.0*ph6_5 - 200.0*ph7_5 + 100.0*ph7_4 + 100.0*ph7_6)^2 + (-200.0*om8_5 + 
	100.0*om9_5 + 100.0*om7_5 - 200.0*om8_5 + 100.0*om8_4 + 100.0*om8_6 + 
	1.5308084989341916e-13*ph9_5 - 1.5308084989341916e-13*ph7_5 - 2500.0*ph8_4 + 
	2500.0*ph8_6)^2 + (-200.0*ps8_5 + 100.0*ps9_5 + 100.0*ps7_5 - 200.0*ps8_5 
	+ 100.0*ps8_4 + 100.0*ps8_6 + 0.25*om8_5)^2 + (-25.0*(ps8_4-ps8_6) * 
	(ph9_5-ph7_5) + 25.0*(ps9_5-ps7_5) * (ph8_4-ph8_6) - 200.0*ph8_5 + 
	100.0*ph9_5 + 100.0*ph7_5 - 200.0*ph8_5 + 100.0*ph8_4 + 100.0*ph8_6)^2 + 
	(-200.0*om9_5 + 100.0*om10_5 + 100.0*om8_5 - 200.0*om9_5 + 100.0*om9_4 + 
	100.0*om9_6 + 1.5308084989341916e-13*ph10_5 - 1.5308084989341916e-13*ph8_5 - 
	2500.0*ph9_4 + 2500.0*ph9_6)^2 + (-200.0*ps9_5 + 100.0*ps10_5 + 
	100.0*ps8_5 - 200.0*ps9_5 + 100.0*ps9_4 + 100.0*ps9_6 + 0.25*om9_5)^2 + 
	(-25.0*(ps9_4-ps9_6) * (ph10_5-ph8_5) + 25.0*(ps10_5-ps8_5) * 
	(ph9_4-ph9_6) - 200.0*ph9_5 + 100.0*ph10_5 + 100.0*ph8_5 - 200.0*ph9_5 + 
	100.0*ph9_4 + 100.0*ph9_6)^2 + (-200.0*omm9_4 + 100.0*omm8_4 + 
	100.0*omm10_4 - 200.0*omm9_4 + 100.0*omm9_3 + 100.0*omm9_5 + 
	1.5308084989341916e-13*phm8_4 - 1.5308084989341916e-13*phm10_4 - 
	2500.0*phm9_3 + 2500.0*phm9_5)^2 + (-200.0*psm9_4 + 100.0*psm8_4 + 
	100.0*psm10_4 - 200.0*psm9_4 + 100.0*psm9_3 + 100.0*psm9_5 + 
	0.25*omm9_4)^2 + (-25.0*(psm9_3-psm9_5) * (phm8_4-phm10_4) + 
	25.0*(psm8_4-psm10_4) * (phm9_3-phm9_5) - 200.0*phm9_4 + 100.0*phm8_4 + 
	100.0*phm10_4 - 200.0*phm9_4 + 100.0*phm9_3 + 100.0*phm9_5)^2 + 
	(-200.0*omm8_4 + 100.0*omm7_4 + 100.0*omm9_4 - 200.0*omm8_4 + 100.0*omm8_3 
	+ 100.0*omm8_5 + 1.5308084989341916e-13*phm7_4 - 
	1.5308084989341916e-13*phm9_4 - 2500.0*phm8_3 + 2500.0*phm8_5)^2 + 
	(-200.0*psm8_4 + 100.0*psm7_4 + 100.0*psm9_4 - 200.0*psm8_4 + 100.0*psm8_3 
	+ 100.0*psm8_5 + 0.25*omm8_4)^2 + (-25.0*(psm8_3-psm8_5) * 
	(phm7_4-phm9_4) + 25.0*(psm7_4-psm9_4) * (phm8_3-phm8_5) - 200.0*phm8_4 
	+ 100.0*phm7_4 + 100.0*phm9_4 - 200.0*phm8_4 + 100.0*phm8_3 + 
	100.0*phm8_5)^2 + (-200.0*omm7_4 + 100.0*omm6_4 + 100.0*omm8_4 - 
	200.0*omm7_4 + 100.0*omm7_3 + 100.0*omm7_5 + 1.5308084989341916e-13*phm6_4 
	- 1.5308084989341916e-13*phm8_4 - 2500.0*phm7_3 + 2500.0*phm7_5)^2 + 
	(-200.0*psm7_4 + 100.0*psm6_4 + 100.0*psm8_4 - 200.0*psm7_4 + 100.0*psm7_3 
	+ 100.0*psm7_5 + 0.25*omm7_4)^2 + (-25.0*(psm7_3-psm7_5) * 
	(phm6_4-phm8_4) + 25.0*(psm6_4-psm8_4) * (phm7_3-phm7_5) - 200.0*phm7_4 
	+ 100.0*phm6_4 + 100.0*phm8_4 - 200.0*phm7_4 + 100.0*phm7_3 + 
	100.0*phm7_5)^2 + (-200.0*omm6_4 + 100.0*omm5_4 + 100.0*omm7_4 - 
	200.0*omm6_4 + 100.0*omm6_3 + 100.0*omm6_5 + 1.5308084989341916e-13*phm5_4 
	- 1.5308084989341916e-13*phm7_4 - 2500.0*phm6_3 + 2500.0*phm6_5)^2 + 
	(-200.0*psm6_4 + 100.0*psm5_4 + 100.0*psm7_4 - 200.0*psm6_4 + 100.0*psm6_3 
	+ 100.0*psm6_5 + 0.25*omm6_4)^2 + (-25.0*(psm6_3-psm6_5) * 
	(phm5_4-phm7_4) + 25.0*(psm5_4-psm7_4) * (phm6_3-phm6_5) - 200.0*phm6_4 
	+ 100.0*phm5_4 + 100.0*phm7_4 - 200.0*phm6_4 + 100.0*phm6_3 + 
	100.0*phm6_5)^2 + (-200.0*omm5_4 + 100.0*omm4_4 + 100.0*omm6_4 - 
	200.0*omm5_4 + 100.0*omm5_3 + 100.0*omm5_5 + 1.5308084989341916e-13*phm4_4 
	- 1.5308084989341916e-13*phm6_4 - 2500.0*phm5_3 + 2500.0*phm5_5)^2 + 
	(-200.0*psm5_4 + 100.0*psm4_4 + 100.0*psm6_4 - 200.0*psm5_4 + 100.0*psm5_3 
	+ 100.0*psm5_5 + 0.25*omm5_4)^2 + (-25.0*(psm5_3-psm5_5) * 
	(phm4_4-phm6_4) + 25.0*(psm4_4-psm6_4) * (phm5_3-phm5_5) - 200.0*phm5_4 
	+ 100.0*phm4_4 + 100.0*phm6_4 - 200.0*phm5_4 + 100.0*phm5_3 + 
	100.0*phm5_5)^2 + (-200.0*omm4_4 + 100.0*omm3_4 + 100.0*omm5_4 - 
	200.0*omm4_4 + 100.0*omm4_3 + 100.0*omm4_5 + 1.5308084989341916e-13*phm3_4 
	- 1.5308084989341916e-13*phm5_4 - 2500.0*phm4_3 + 2500.0*phm4_5)^2 + 
	(-200.0*psm4_4 + 100.0*psm3_4 + 100.0*psm5_4 - 200.0*psm4_4 + 100.0*psm4_3 
	+ 100.0*psm4_5 + 0.25*omm4_4)^2 + (-25.0*(psm4_3-psm4_5) * 
	(phm3_4-phm5_4) + 25.0*(psm3_4-psm5_4) * (phm4_3-phm4_5) - 200.0*phm4_4 
	+ 100.0*phm3_4 + 100.0*phm5_4 - 200.0*phm4_4 + 100.0*phm4_3 + 
	100.0*phm4_5)^2 + (-200.0*omm3_4 + 100.0*omm2_4 + 100.0*omm4_4 - 
	200.0*omm3_4 + 100.0*omm3_3 + 100.0*omm3_5 + 1.5308084989341916e-13*phm2_4 
	- 1.5308084989341916e-13*phm4_4 - 2500.0*phm3_3 + 2500.0*phm3_5)^2 + 
	(-200.0*psm3_4 + 100.0*psm2_4 + 100.0*psm4_4 - 200.0*psm3_4 + 100.0*psm3_3 
	+ 100.0*psm3_5 + 0.25*omm3_4)^2 + (-25.0*(psm3_3-psm3_5) * 
	(phm2_4-phm4_4) + 25.0*(psm2_4-psm4_4) * (phm3_3-phm3_5) - 200.0*phm3_4 
	+ 100.0*phm2_4 + 100.0*phm4_4 - 200.0*phm3_4 + 100.0*phm3_3 + 
	100.0*phm3_5)^2 + (-200.0*omm2_4 + 100.0*omm1_4 + 100.0*omm3_4 - 
	200.0*omm2_4 + 100.0*omm2_3 + 100.0*omm2_5 + 1.5308084989341916e-13*phm1_4 
	- 1.5308084989341916e-13*phm3_4 - 2500.0*phm2_3 + 2500.0*phm2_5)^2 + 
	(-200.0*psm2_4 + 100.0*psm1_4 + 100.0*psm3_4 - 200.0*psm2_4 + 100.0*psm2_3 
	+ 100.0*psm2_5 + 0.25*omm2_4)^2 + (-25.0*(psm2_3-psm2_5) * 
	(phm1_4-phm3_4) + 25.0*(psm1_4-psm3_4) * (phm2_3-phm2_5) - 200.0*phm2_4 
	+ 100.0*phm1_4 + 100.0*phm3_4 - 200.0*phm2_4 + 100.0*phm2_3 + 
	100.0*phm2_5)^2 + (-200.0*omm1_4 + 100.0*om0_4 + 100.0*omm2_4 - 
	200.0*omm1_4 + 100.0*omm1_3 + 100.0*omm1_5 + 1.5308084989341916e-13*ph0_4 - 
	1.5308084989341916e-13*phm2_4 - 2500.0*phm1_3 + 2500.0*phm1_5)^2 + 
	(-200.0*psm1_4 + 100.0*ps0_4 + 100.0*psm2_4 - 200.0*psm1_4 + 100.0*psm1_3 
	+ 100.0*psm1_5 + 0.25*omm1_4)^2 + (-25.0*(psm1_3-psm1_5) * (ph0_4-phm2_4) 
	+ 25.0*(ps0_4-psm2_4) * (phm1_3-phm1_5) - 200.0*phm1_4 + 100.0*ph0_4 + 
	100.0*phm2_4 - 200.0*phm1_4 + 100.0*phm1_3 + 100.0*phm1_5)^2 + 
	(-200.0*om0_4 + 100.0*om1_4 + 100.0*omm1_4 - 200.0*om0_4 + 100.0*om0_3 + 
	100.0*om0_5 + 1.5308084989341916e-13*ph1_4 - 1.5308084989341916e-13*phm1_4 - 
	2500.0*ph0_3 + 2500.0*ph0_5)^2 + (-200.0*ps0_4 + 100.0*ps1_4 + 
	100.0*psm1_4 - 200.0*ps0_4 + 100.0*ps0_3 + 100.0*ps0_5 + 0.25*om0_4)^2 + 
	(-25.0*(ps0_3-ps0_5) * (ph1_4-phm1_4) + 25.0*(ps1_4-psm1_4) * 
	(ph0_3-ph0_5) - 200.0*ph0_4 + 100.0*ph1_4 + 100.0*phm1_4 - 200.0*ph0_4 + 
	100.0*ph0_3 + 100.0*ph0_5)^2 + (-200.0*om1_4 + 100.0*om2_4 + 100.0*om0_4 - 
	200.0*om1_4 + 100.0*om1_3 + 100.0*om1_5 + 1.5308084989341916e-13*ph2_4 - 
	1.5308084989341916e-13*ph0_4 - 2500.0*ph1_3 + 2500.0*ph1_5)^2 + 
	(-200.0*ps1_4 + 100.0*ps2_4 + 100.0*ps0_4 - 200.0*ps1_4 + 100.0*ps1_3 + 
	100.0*ps1_5 + 0.25*om1_4)^2 + (-25.0*(ps1_3-ps1_5) * (ph2_4-ph0_4) + 
	25.0*(ps2_4-ps0_4) * (ph1_3-ph1_5) - 200.0*ph1_4 + 100.0*ph2_4 + 
	100.0*ph0_4 - 200.0*ph1_4 + 100.0*ph1_3 + 100.0*ph1_5)^2 + (-200.0*om2_4 + 
	100.0*om3_4 + 100.0*om1_4 - 200.0*om2_4 + 100.0*om2_3 + 100.0*om2_5 + 
	1.5308084989341916e-13*ph3_4 - 1.5308084989341916e-13*ph1_4 - 2500.0*ph2_3 + 
	2500.0*ph2_5)^2 + (-200.0*ps2_4 + 100.0*ps3_4 + 100.0*ps1_4 - 200.0*ps2_4 
	+ 100.0*ps2_3 + 100.0*ps2_5 + 0.25*om2_4)^2 + (-25.0*(ps2_3-ps2_5) * 
	(ph3_4-ph1_4) + 25.0*(ps3_4-ps1_4) * (ph2_3-ph2_5) - 200.0*ph2_4 + 
	100.0*ph3_4 + 100.0*ph1_4 - 200.0*ph2_4 + 100.0*ph2_3 + 100.0*ph2_5)^2 + 
	(-200.0*om3_4 + 100.0*om4_4 + 100.0*om2_4 - 200.0*om3_4 + 100.0*om3_3 + 
	100.0*om3_5 + 1.5308084989341916e-13*ph4_4 - 1.5308084989341916e-13*ph2_4 - 
	2500.0*ph3_3 + 2500.0*ph3_5)^2 + (-200.0*ps3_4 + 100.0*ps4_4 + 100.0*ps2_4 
	- 200.0*ps3_4 + 100.0*ps3_3 + 100.0*ps3_5 + 0.25*om3_4)^2 + 
	(-25.0*(ps3_3-ps3_5) * (ph4_4-ph2_4) + 25.0*(ps4_4-ps2_4) * 
	(ph3_3-ph3_5) - 200.0*ph3_4 + 100.0*ph4_4 + 100.0*ph2_4 - 200.0*ph3_4 + 
	100.0*ph3_3 + 100.0*ph3_5)^2 + (-200.0*om4_4 + 100.0*om5_4 + 100.0*om3_4 - 
	200.0*om4_4 + 100.0*om4_3 + 100.0*om4_5 + 1.5308084989341916e-13*ph5_4 - 
	1.5308084989341916e-13*ph3_4 - 2500.0*ph4_3 + 2500.0*ph4_5)^2 + 
	(-200.0*ps4_4 + 100.0*ps5_4 + 100.0*ps3_4 - 200.0*ps4_4 + 100.0*ps4_3 + 
	100.0*ps4_5 + 0.25*om4_4)^2 + (-25.0*(ps4_3-ps4_5) * (ph5_4-ph3_4) + 
	25.0*(ps5_4-ps3_4) * (ph4_3-ph4_5) - 200.0*ph4_4 + 100.0*ph5_4 + 
	100.0*ph3_4 - 200.0*ph4_4 + 100.0*ph4_3 + 100.0*ph4_5)^2 + (-200.0*om5_4 + 
	100.0*om6_4 + 100.0*om4_4 - 200.0*om5_4 + 100.0*om5_3 + 100.0*om5_5 + 
	1.5308084989341916e-13*ph6_4 - 1.5308084989341916e-13*ph4_4 - 2500.0*ph5_3 + 
	2500.0*ph5_5)^2 + (-200.0*ps5_4 + 100.0*ps6_4 + 100.0*ps4_4 - 200.0*ps5_4 
	+ 100.0*ps5_3 + 100.0*ps5_5 + 0.25*om5_4)^2 + (-25.0*(ps5_3-ps5_5) * 
	(ph6_4-ph4_4) + 25.0*(ps6_4-ps4_4) * (ph5_3-ph5_5) - 200.0*ph5_4 + 
	100.0*ph6_4 + 100.0*ph4_4 - 200.0*ph5_4 + 100.0*ph5_3 + 100.0*ph5_5)^2 + 
	(-200.0*om6_4 + 100.0*om7_4 + 100.0*om5_4 - 200.0*om6_4 + 100.0*om6_3 + 
	100.0*om6_5 + 1.5308084989341916e-13*ph7_4 - 1.5308084989341916e-13*ph5_4 - 
	2500.0*ph6_3 + 2500.0*ph6_5)^2 + (-200.0*ps6_4 + 100.0*ps7_4 + 100.0*ps5_4 
	- 200.0*ps6_4 + 100.0*ps6_3 + 100.0*ps6_5 + 0.25*om6_4)^2 + 
	(-25.0*(ps6_3-ps6_5) * (ph7_4-ph5_4) + 25.0*(ps7_4-ps5_4) * 
	(ph6_3-ph6_5) - 200.0*ph6_4 + 100.0*ph7_4 + 100.0*ph5_4 - 200.0*ph6_4 + 
	100.0*ph6_3 + 100.0*ph6_5)^2 + (-200.0*om7_4 + 100.0*om8_4 + 100.0*om6_4 - 
	200.0*om7_4 + 100.0*om7_3 + 100.0*om7_5 + 1.5308084989341916e-13*ph8_4 - 
	1.5308084989341916e-13*ph6_4 - 2500.0*ph7_3 + 2500.0*ph7_5)^2 + 
	(-200.0*ps7_4 + 100.0*ps8_4 + 100.0*ps6_4 - 200.0*ps7_4 + 100.0*ps7_3 + 
	100.0*ps7_5 + 0.25*om7_4)^2 + (-25.0*(ps7_3-ps7_5) * (ph8_4-ph6_4) + 
	25.0*(ps8_4-ps6_4) * (ph7_3-ph7_5) - 200.0*ph7_4 + 100.0*ph8_4 + 
	100.0*ph6_4 - 200.0*ph7_4 + 100.0*ph7_3 + 100.0*ph7_5)^2 + (-200.0*om8_4 + 
	100.0*om9_4 + 100.0*om7_4 - 200.0*om8_4 + 100.0*om8_3 + 100.0*om8_5 + 
	1.5308084989341916e-13*ph9_4 - 1.5308084989341916e-13*ph7_4 - 2500.0*ph8_3 + 
	2500.0*ph8_5)^2 + (-200.0*ps8_4 + 100.0*ps9_4 + 100.0*ps7_4 - 200.0*ps8_4 
	+ 100.0*ps8_3 + 100.0*ps8_5 + 0.25*om8_4)^2 + (-25.0*(ps8_3-ps8_5) * 
	(ph9_4-ph7_4) + 25.0*(ps9_4-ps7_4) * (ph8_3-ph8_5) - 200.0*ph8_4 + 
	100.0*ph9_4 + 100.0*ph7_4 - 200.0*ph8_4 + 100.0*ph8_3 + 100.0*ph8_5)^2 + 
	(-200.0*om9_4 + 100.0*om10_4 + 100.0*om8_4 - 200.0*om9_4 + 100.0*om9_3 + 
	100.0*om9_5 + 1.5308084989341916e-13*ph10_4 - 1.5308084989341916e-13*ph8_4 - 
	2500.0*ph9_3 + 2500.0*ph9_5)^2 + (-200.0*ps9_4 + 100.0*ps10_4 + 
	100.0*ps8_4 - 200.0*ps9_4 + 100.0*ps9_3 + 100.0*ps9_5 + 0.25*om9_4)^2 + 
	(-25.0*(ps9_3-ps9_5) * (ph10_4-ph8_4) + 25.0*(ps10_4-ps8_4) * 
	(ph9_3-ph9_5) - 200.0*ph9_4 + 100.0*ph10_4 + 100.0*ph8_4 - 200.0*ph9_4 + 
	100.0*ph9_3 + 100.0*ph9_5)^2 + (-200.0*omm9_3 + 100.0*omm8_3 + 
	100.0*omm10_3 - 200.0*omm9_3 + 100.0*omm9_2 + 100.0*omm9_4 + 
	1.5308084989341916e-13*phm8_3 - 1.5308084989341916e-13*phm10_3 - 
	2500.0*phm9_2 + 2500.0*phm9_4)^2 + (-200.0*psm9_3 + 100.0*psm8_3 + 
	100.0*psm10_3 - 200.0*psm9_3 + 100.0*psm9_2 + 100.0*psm9_4 + 
	0.25*omm9_3)^2 + (-25.0*(psm9_2-psm9_4) * (phm8_3-phm10_3) + 
	25.0*(psm8_3-psm10_3) * (phm9_2-phm9_4) - 200.0*phm9_3 + 100.0*phm8_3 + 
	100.0*phm10_3 - 200.0*phm9_3 + 100.0*phm9_2 + 100.0*phm9_4)^2 + 
	(-200.0*omm8_3 + 100.0*omm7_3 + 100.0*omm9_3 - 200.0*omm8_3 + 100.0*omm8_2 
	+ 100.0*omm8_4 + 1.5308084989341916e-13*phm7_3 - 
	1.5308084989341916e-13*phm9_3 - 2500.0*phm8_2 + 2500.0*phm8_4)^2 + 
	(-200.0*psm8_3 + 100.0*psm7_3 + 100.0*psm9_3 - 200.0*psm8_3 + 100.0*psm8_2 
	+ 100.0*psm8_4 + 0.25*omm8_3)^2 + (-25.0*(psm8_2-psm8_4) * 
	(phm7_3-phm9_3) + 25.0*(psm7_3-psm9_3) * (phm8_2-phm8_4) - 200.0*phm8_3 
	+ 100.0*phm7_3 + 100.0*phm9_3 - 200.0*phm8_3 + 100.0*phm8_2 + 
	100.0*phm8_4)^2 + (-200.0*omm7_3 + 100.0*omm6_3 + 100.0*omm8_3 - 
	200.0*omm7_3 + 100.0*omm7_2 + 100.0*omm7_4 + 1.5308084989341916e-13*phm6_3 
	- 1.5308084989341916e-13*phm8_3 - 2500.0*phm7_2 + 2500.0*phm7_4)^2 + 
	(-200.0*psm7_3 + 100.0*psm6_3 + 100.0*psm8_3 - 200.0*psm7_3 + 100.0*psm7_2 
	+ 100.0*psm7_4 + 0.25*omm7_3)^2 + (-25.0*(psm7_2-psm7_4) * 
	(phm6_3-phm8_3) + 25.0*(psm6_3-psm8_3) * (phm7_2-phm7_4) - 200.0*phm7_3 
	+ 100.0*phm6_3 + 100.0*phm8_3 - 200.0*phm7_3 + 100.0*phm7_2 + 
	100.0*phm7_4)^2 + (-200.0*omm6_3 + 100.0*omm5_3 + 100.0*omm7_3 - 
	200.0*omm6_3 + 100.0*omm6_2 + 100.0*omm6_4 + 1.5308084989341916e-13*phm5_3 
	- 1.5308084989341916e-13*phm7_3 - 2500.0*phm6_2 + 2500.0*phm6_4)^2 + 
	(-200.0*psm6_3 + 100.0*psm5_3 + 100.0*psm7_3 - 200.0*psm6_3 + 100.0*psm6_2 
	+ 100.0*psm6_4 + 0.25*omm6_3)^2 + (-25.0*(psm6_2-psm6_4) * 
	(phm5_3-phm7_3) + 25.0*(psm5_3-psm7_3) * (phm6_2-phm6_4) - 200.0*phm6_3 
	+ 100.0*phm5_3 + 100.0*phm7_3 - 200.0*phm6_3 + 100.0*phm6_2 + 
	100.0*phm6_4)^2 + (-200.0*omm5_3 + 100.0*omm4_3 + 100.0*omm6_3 - 
	200.0*omm5_3 + 100.0*omm5_2 + 100.0*omm5_4 + 1.5308084989341916e-13*phm4_3 
	- 1.5308084989341916e-13*phm6_3 - 2500.0*phm5_2 + 2500.0*phm5_4)^2 + 
	(-200.0*psm5_3 + 100.0*psm4_3 + 100.0*psm6_3 - 200.0*psm5_3 + 100.0*psm5_2 
	+ 100.0*psm5_4 + 0.25*omm5_3)^2 + (-25.0*(psm5_2-psm5_4) * 
	(phm4_3-phm6_3) + 25.0*(psm4_3-psm6_3) * (phm5_2-phm5_4) - 200.0*phm5_3 
	+ 100.0*phm4_3 + 100.0*phm6_3 - 200.0*phm5_3 + 100.0*phm5_2 + 
	100.0*phm5_4)^2 + (-200.0*omm4_3 + 100.0*omm3_3 + 100.0*omm5_3 - 
	200.0*omm4_3 + 100.0*omm4_2 + 100.0*omm4_4 + 1.5308084989341916e-13*phm3_3 
	- 1.5308084989341916e-13*phm5_3 - 2500.0*phm4_2 + 2500.0*phm4_4)^2 + 
	(-200.0*psm4_3 + 100.0*psm3_3 + 100.0*psm5_3 - 200.0*psm4_3 + 100.0*psm4_2 
	+ 100.0*psm4_4 + 0.25*omm4_3)^2 + (-25.0*(psm4_2-psm4_4) * 
	(phm3_3-phm5_3) + 25.0*(psm3_3-psm5_3) * (phm4_2-phm4_4) - 200.0*phm4_3 
	+ 100.0*phm3_3 + 100.0*phm5_3 - 200.0*phm4_3 + 100.0*phm4_2 + 
	100.0*phm4_4)^2 + (-200.0*omm3_3 + 100.0*omm2_3 + 100.0*omm4_3 - 
	200.0*omm3_3 + 100.0*omm3_2 + 100.0*omm3_4 + 1.5308084989341916e-13*phm2_3 
	- 1.5308084989341916e-13*phm4_3 - 2500.0*phm3_2 + 2500.0*phm3_4)^2 + 
	(-200.0*psm3_3 + 100.0*psm2_3 + 100.0*psm4_3 - 200.0*psm3_3 + 100.0*psm3_2 
	+ 100.0*psm3_4 + 0.25*omm3_3)^2 + (-25.0*(psm3_2-psm3_4) * 
	(phm2_3-phm4_3) + 25.0*(psm2_3-psm4_3) * (phm3_2-phm3_4) - 200.0*phm3_3 
	+ 100.0*phm2_3 + 100.0*phm4_3 - 200.0*phm3_3 + 100.0*phm3_2 + 
	100.0*phm3_4)^2 + (-200.0*omm2_3 + 100.0*omm1_3 + 100.0*omm3_3 - 
	200.0*omm2_3 + 100.0*omm2_2 + 100.0*omm2_4 + 1.5308084989341916e-13*phm1_3 
	- 1.5308084989341916e-13*phm3_3 - 2500.0*phm2_2 + 2500.0*phm2_4)^2 + 
	(-200.0*psm2_3 + 100.0*psm1_3 + 100.0*psm3_3 - 200.0*psm2_3 + 100.0*psm2_2 
	+ 100.0*psm2_4 + 0.25*omm2_3)^2 + (-25.0*(psm2_2-psm2_4) * 
	(phm1_3-phm3_3) + 25.0*(psm1_3-psm3_3) * (phm2_2-phm2_4) - 200.0*phm2_3 
	+ 100.0*phm1_3 + 100.0*phm3_3 - 200.0*phm2_3 + 100.0*phm2_2 + 
	100.0*phm2_4)^2 + (-200.0*omm1_3 + 100.0*om0_3 + 100.0*omm2_3 - 
	200.0*omm1_3 + 100.0*omm1_2 + 100.0*omm1_4 + 1.5308084989341916e-13*ph0_3 - 
	1.5308084989341916e-13*phm2_3 - 2500.0*phm1_2 + 2500.0*phm1_4)^2 + 
	(-200.0*psm1_3 + 100.0*ps0_3 + 100.0*psm2_3 - 200.0*psm1_3 + 100.0*psm1_2 
	+ 100.0*psm1_4 + 0.25*omm1_3)^2 + (-25.0*(psm1_2-psm1_4) * (ph0_3-phm2_3) 
	+ 25.0*(ps0_3-psm2_3) * (phm1_2-phm1_4) - 200.0*phm1_3 + 100.0*ph0_3 + 
	100.0*phm2_3 - 200.0*phm1_3 + 100.0*phm1_2 + 100.0*phm1_4)^2 + 
	(-200.0*om0_3 + 100.0*om1_3 + 100.0*omm1_3 - 200.0*om0_3 + 100.0*om0_2 + 
	100.0*om0_4 + 1.5308084989341916e-13*ph1_3 - 1.5308084989341916e-13*phm1_3 - 
	2500.0*ph0_2 + 2500.0*ph0_4)^2 + (-200.0*ps0_3 + 100.0*ps1_3 + 
	100.0*psm1_3 - 200.0*ps0_3 + 100.0*ps0_2 + 100.0*ps0_4 + 0.25*om0_3)^2 + 
	(-25.0*(ps0_2-ps0_4) * (ph1_3-phm1_3) + 25.0*(ps1_3-psm1_3) * 
	(ph0_2-ph0_4) - 200.0*ph0_3 + 100.0*ph1_3 + 100.0*phm1_3 - 200.0*ph0_3 + 
	100.0*ph0_2 + 100.0*ph0_4)^2 + (-200.0*om1_3 + 100.0*om2_3 + 100.0*om0_3 - 
	200.0*om1_3 + 100.0*om1_2 + 100.0*om1_4 + 1.5308084989341916e-13*ph2_3 - 
	1.5308084989341916e-13*ph0_3 - 2500.0*ph1_2 + 2500.0*ph1_4)^2 + 
	(-200.0*ps1_3 + 100.0*ps2_3 + 100.0*ps0_3 - 200.0*ps1_3 + 100.0*ps1_2 + 
	100.0*ps1_4 + 0.25*om1_3)^2 + (-25.0*(ps1_2-ps1_4) * (ph2_3-ph0_3) + 
	25.0*(ps2_3-ps0_3) * (ph1_2-ph1_4) - 200.0*ph1_3 + 100.0*ph2_3 + 
	100.0*ph0_3 - 200.0*ph1_3 + 100.0*ph1_2 + 100.0*ph1_4)^2 + (-200.0*om2_3 + 
	100.0*om3_3 + 100.0*om1_3 - 200.0*om2_3 + 100.0*om2_2 + 100.0*om2_4 + 
	1.5308084989341916e-13*ph3_3 - 1.5308084989341916e-13*ph1_3 - 2500.0*ph2_2 + 
	2500.0*ph2_4)^2 + (-200.0*ps2_3 + 100.0*ps3_3 + 100.0*ps1_3 - 200.0*ps2_3 
	+ 100.0*ps2_2 + 100.0*ps2_4 + 0.25*om2_3)^2 + (-25.0*(ps2_2-ps2_4) * 
	(ph3_3-ph1_3) + 25.0*(ps3_3-ps1_3) * (ph2_2-ph2_4) - 200.0*ph2_3 + 
	100.0*ph3_3 + 100.0*ph1_3 - 200.0*ph2_3 + 100.0*ph2_2 + 100.0*ph2_4)^2 + 
	(-200.0*om3_3 + 100.0*om4_3 + 100.0*om2_3 - 200.0*om3_3 + 100.0*om3_2 + 
	100.0*om3_4 + 1.5308084989341916e-13*ph4_3 - 1.5308084989341916e-13*ph2_3 - 
	2500.0*ph3_2 + 2500.0*ph3_4)^2 + (-200.0*ps3_3 + 100.0*ps4_3 + 100.0*ps2_3 
	- 200.0*ps3_3 + 100.0*ps3_2 + 100.0*ps3_4 + 0.25*om3_3)^2 + 
	(-25.0*(ps3_2-ps3_4) * (ph4_3-ph2_3) + 25.0*(ps4_3-ps2_3) * 
	(ph3_2-ph3_4) - 200.0*ph3_3 + 100.0*ph4_3 + 100.0*ph2_3 - 200.0*ph3_3 + 
	100.0*ph3_2 + 100.0*ph3_4)^2 + (-200.0*om4_3 + 100.0*om5_3 + 100.0*om3_3 - 
	200.0*om4_3 + 100.0*om4_2 + 100.0*om4_4 + 1.5308084989341916e-13*ph5_3 - 
	1.5308084989341916e-13*ph3_3 - 2500.0*ph4_2 + 2500.0*ph4_4)^2 + 
	(-200.0*ps4_3 + 100.0*ps5_3 + 100.0*ps3_3 - 200.0*ps4_3 + 100.0*ps4_2 + 
	100.0*ps4_4 + 0.25*om4_3)^2 + (-25.0*(ps4_2-ps4_4) * (ph5_3-ph3_3) + 
	25.0*(ps5_3-ps3_3) * (ph4_2-ph4_4) - 200.0*ph4_3 + 100.0*ph5_3 + 
	100.0*ph3_3 - 200.0*ph4_3 + 100.0*ph4_2 + 100.0*ph4_4)^2 + (-200.0*om5_3 + 
	100.0*om6_3 + 100.0*om4_3 - 200.0*om5_3 + 100.0*om5_2 + 100.0*om5_4 + 
	1.5308084989341916e-13*ph6_3 - 1.5308084989341916e-13*ph4_3 - 2500.0*ph5_2 + 
	2500.0*ph5_4)^2 + (-200.0*ps5_3 + 100.0*ps6_3 + 100.0*ps4_3 - 200.0*ps5_3 
	+ 100.0*ps5_2 + 100.0*ps5_4 + 0.25*om5_3)^2 + (-25.0*(ps5_2-ps5_4) * 
	(ph6_3-ph4_3) + 25.0*(ps6_3-ps4_3) * (ph5_2-ph5_4) - 200.0*ph5_3 + 
	100.0*ph6_3 + 100.0*ph4_3 - 200.0*ph5_3 + 100.0*ph5_2 + 100.0*ph5_4)^2 + 
	(-200.0*om6_3 + 100.0*om7_3 + 100.0*om5_3 - 200.0*om6_3 + 100.0*om6_2 + 
	100.0*om6_4 + 1.5308084989341916e-13*ph7_3 - 1.5308084989341916e-13*ph5_3 - 
	2500.0*ph6_2 + 2500.0*ph6_4)^2 + (-200.0*ps6_3 + 100.0*ps7_3 + 100.0*ps5_3 
	- 200.0*ps6_3 + 100.0*ps6_2 + 100.0*ps6_4 + 0.25*om6_3)^2 + 
	(-25.0*(ps6_2-ps6_4) * (ph7_3-ph5_3) + 25.0*(ps7_3-ps5_3) * 
	(ph6_2-ph6_4) - 200.0*ph6_3 + 100.0*ph7_3 + 100.0*ph5_3 - 200.0*ph6_3 + 
	100.0*ph6_2 + 100.0*ph6_4)^2 + (-200.0*om7_3 + 100.0*om8_3 + 100.0*om6_3 - 
	200.0*om7_3 + 100.0*om7_2 + 100.0*om7_4 + 1.5308084989341916e-13*ph8_3 - 
	1.5308084989341916e-13*ph6_3 - 2500.0*ph7_2 + 2500.0*ph7_4)^2 + 
	(-200.0*ps7_3 + 100.0*ps8_3 + 100.0*ps6_3 - 200.0*ps7_3 + 100.0*ps7_2 + 
	100.0*ps7_4 + 0.25*om7_3)^2 + (-25.0*(ps7_2-ps7_4) * (ph8_3-ph6_3) + 
	25.0*(ps8_3-ps6_3) * (ph7_2-ph7_4) - 200.0*ph7_3 + 100.0*ph8_3 + 
	100.0*ph6_3 - 200.0*ph7_3 + 100.0*ph7_2 + 100.0*ph7_4)^2 + (-200.0*om8_3 + 
	100.0*om9_3 + 100.0*om7_3 - 200.0*om8_3 + 100.0*om8_2 + 100.0*om8_4 + 
	1.5308084989341916e-13*ph9_3 - 1.5308084989341916e-13*ph7_3 - 2500.0*ph8_2 + 
	2500.0*ph8_4)^2 + (-200.0*ps8_3 + 100.0*ps9_3 + 100.0*ps7_3 - 200.0*ps8_3 
	+ 100.0*ps8_2 + 100.0*ps8_4 + 0.25*om8_3)^2 + (-25.0*(ps8_2-ps8_4) * 
	(ph9_3-ph7_3) + 25.0*(ps9_3-ps7_3) * (ph8_2-ph8_4) - 200.0*ph8_3 + 
	100.0*ph9_3 + 100.0*ph7_3 - 200.0*ph8_3 + 100.0*ph8_2 + 100.0*ph8_4)^2 + 
	(-200.0*om9_3 + 100.0*om10_3 + 100.0*om8_3 - 200.0*om9_3 + 100.0*om9_2 + 
	100.0*om9_4 + 1.5308084989341916e-13*ph10_3 - 1.5308084989341916e-13*ph8_3 - 
	2500.0*ph9_2 + 2500.0*ph9_4)^2 + (-200.0*ps9_3 + 100.0*ps10_3 + 
	100.0*ps8_3 - 200.0*ps9_3 + 100.0*ps9_2 + 100.0*ps9_4 + 0.25*om9_3)^2 + 
	(-25.0*(ps9_2-ps9_4) * (ph10_3-ph8_3) + 25.0*(ps10_3-ps8_3) * 
	(ph9_2-ph9_4) - 200.0*ph9_3 + 100.0*ph10_3 + 100.0*ph8_3 - 200.0*ph9_3 + 
	100.0*ph9_2 + 100.0*ph9_4)^2 + (-200.0*omm9_2 + 100.0*omm8_2 + 
	100.0*omm10_2 - 200.0*omm9_2 + 100.0*omm9_1 + 100.0*omm9_3 + 
	1.5308084989341916e-13*phm8_2 - 1.5308084989341916e-13*phm10_2 - 
	2500.0*phm9_1 + 2500.0*phm9_3)^2 + (-200.0*psm9_2 + 100.0*psm8_2 + 
	100.0*psm10_2 - 200.0*psm9_2 + 100.0*psm9_1 + 100.0*psm9_3 + 
	0.25*omm9_2)^2 + (-25.0*(psm9_1-psm9_3) * (phm8_2-phm10_2) + 
	25.0*(psm8_2-psm10_2) * (phm9_1-phm9_3) - 200.0*phm9_2 + 100.0*phm8_2 + 
	100.0*phm10_2 - 200.0*phm9_2 + 100.0*phm9_1 + 100.0*phm9_3)^2 + 
	(-200.0*omm8_2 + 100.0*omm7_2 + 100.0*omm9_2 - 200.0*omm8_2 + 100.0*omm8_1 
	+ 100.0*omm8_3 + 1.5308084989341916e-13*phm7_2 - 
	1.5308084989341916e-13*phm9_2 - 2500.0*phm8_1 + 2500.0*phm8_3)^2 + 
	(-200.0*psm8_2 + 100.0*psm7_2 + 100.0*psm9_2 - 200.0*psm8_2 + 100.0*psm8_1 
	+ 100.0*psm8_3 + 0.25*omm8_2)^2 + (-25.0*(psm8_1-psm8_3) * 
	(phm7_2-phm9_2) + 25.0*(psm7_2-psm9_2) * (phm8_1-phm8_3) - 200.0*phm8_2 
	+ 100.0*phm7_2 + 100.0*phm9_2 - 200.0*phm8_2 + 100.0*phm8_1 + 
	100.0*phm8_3)^2 + (-200.0*omm7_2 + 100.0*omm6_2 + 100.0*omm8_2 - 
	200.0*omm7_2 + 100.0*omm7_1 + 100.0*omm7_3 + 1.5308084989341916e-13*phm6_2 
	- 1.5308084989341916e-13*phm8_2 - 2500.0*phm7_1 + 2500.0*phm7_3)^2 + 
	(-200.0*psm7_2 + 100.0*psm6_2 + 100.0*psm8_2 - 200.0*psm7_2 + 100.0*psm7_1 
	+ 100.0*psm7_3 + 0.25*omm7_2)^2 + (-25.0*(psm7_1-psm7_3) * 
	(phm6_2-phm8_2) + 25.0*(psm6_2-psm8_2) * (phm7_1-phm7_3) - 200.0*phm7_2 
	+ 100.0*phm6_2 + 100.0*phm8_2 - 200.0*phm7_2 + 100.0*phm7_1 + 
	100.0*phm7_3)^2 + (-200.0*omm6_2 + 100.0*omm5_2 + 100.0*omm7_2 - 
	200.0*omm6_2 + 100.0*omm6_1 + 100.0*omm6_3 + 1.5308084989341916e-13*phm5_2 
	- 1.5308084989341916e-13*phm7_2 - 2500.0*phm6_1 + 2500.0*phm6_3)^2 + 
	(-200.0*psm6_2 + 100.0*psm5_2 + 100.0*psm7_2 - 200.0*psm6_2 + 100.0*psm6_1 
	+ 100.0*psm6_3 + 0.25*omm6_2)^2 + (-25.0*(psm6_1-psm6_3) * 
	(phm5_2-phm7_2) + 25.0*(psm5_2-psm7_2) * (phm6_1-phm6_3) - 200.0*phm6_2 
	+ 100.0*phm5_2 + 100.0*phm7_2 - 200.0*phm6_2 + 100.0*phm6_1 + 
	100.0*phm6_3)^2 + (-200.0*omm5_2 + 100.0*omm4_2 + 100.0*omm6_2 - 
	200.0*omm5_2 + 100.0*omm5_1 + 100.0*omm5_3 + 1.5308084989341916e-13*phm4_2 
	- 1.5308084989341916e-13*phm6_2 - 2500.0*phm5_1 + 2500.0*phm5_3)^2 + 
	(-200.0*psm5_2 + 100.0*psm4_2 + 100.0*psm6_2 - 200.0*psm5_2 + 100.0*psm5_1 
	+ 100.0*psm5_3 + 0.25*omm5_2)^2 + (-25.0*(psm5_1-psm5_3) * 
	(phm4_2-phm6_2) + 25.0*(psm4_2-psm6_2) * (phm5_1-phm5_3) - 200.0*phm5_2 
	+ 100.0*phm4_2 + 100.0*phm6_2 - 200.0*phm5_2 + 100.0*phm5_1 + 
	100.0*phm5_3)^2 + (-200.0*omm4_2 + 100.0*omm3_2 + 100.0*omm5_2 - 
	200.0*omm4_2 + 100.0*omm4_1 + 100.0*omm4_3 + 1.5308084989341916e-13*phm3_2 
	- 1.5308084989341916e-13*phm5_2 - 2500.0*phm4_1 + 2500.0*phm4_3)^2 + 
	(-200.0*psm4_2 + 100.0*psm3_2 + 100.0*psm5_2 - 200.0*psm4_2 + 100.0*psm4_1 
	+ 100.0*psm4_3 + 0.25*omm4_2)^2 + (-25.0*(psm4_1-psm4_3) * 
	(phm3_2-phm5_2) + 25.0*(psm3_2-psm5_2) * (phm4_1-phm4_3) - 200.0*phm4_2 
	+ 100.0*phm3_2 + 100.0*phm5_2 - 200.0*phm4_2 + 100.0*phm4_1 + 
	100.0*phm4_3)^2 + (-200.0*omm3_2 + 100.0*omm2_2 + 100.0*omm4_2 - 
	200.0*omm3_2 + 100.0*omm3_1 + 100.0*omm3_3 + 1.5308084989341916e-13*phm2_2 
	- 1.5308084989341916e-13*phm4_2 - 2500.0*phm3_1 + 2500.0*phm3_3)^2 + 
	(-200.0*psm3_2 + 100.0*psm2_2 + 100.0*psm4_2 - 200.0*psm3_2 + 100.0*psm3_1 
	+ 100.0*psm3_3 + 0.25*omm3_2)^2 + (-25.0*(psm3_1-psm3_3) * 
	(phm2_2-phm4_2) + 25.0*(psm2_2-psm4_2) * (phm3_1-phm3_3) - 200.0*phm3_2 
	+ 100.0*phm2_2 + 100.0*phm4_2 - 200.0*phm3_2 + 100.0*phm3_1 + 
	100.0*phm3_3)^2 + (-200.0*omm2_2 + 100.0*omm1_2 + 100.0*omm3_2 - 
	200.0*omm2_2 + 100.0*omm2_1 + 100.0*omm2_3 + 1.5308084989341916e-13*phm1_2 
	- 1.5308084989341916e-13*phm3_2 - 2500.0*phm2_1 + 2500.0*phm2_3)^2 + 
	(-200.0*psm2_2 + 100.0*psm1_2 + 100.0*psm3_2 - 200.0*psm2_2 + 100.0*psm2_1 
	+ 100.0*psm2_3 + 0.25*omm2_2)^2 + (-25.0*(psm2_1-psm2_3) * 
	(phm1_2-phm3_2) + 25.0*(psm1_2-psm3_2) * (phm2_1-phm2_3) - 200.0*phm2_2 
	+ 100.0*phm1_2 + 100.0*phm3_2 - 200.0*phm2_2 + 100.0*phm2_1 + 
	100.0*phm2_3)^2 + (-200.0*omm1_2 + 100.0*om0_2 + 100.0*omm2_2 - 
	200.0*omm1_2 + 100.0*omm1_1 + 100.0*omm1_3 + 1.5308084989341916e-13*ph0_2 - 
	1.5308084989341916e-13*phm2_2 - 2500.0*phm1_1 + 2500.0*phm1_3)^2 + 
	(-200.0*psm1_2 + 100.0*ps0_2 + 100.0*psm2_2 - 200.0*psm1_2 + 100.0*psm1_1 
	+ 100.0*psm1_3 + 0.25*omm1_2)^2 + (-25.0*(psm1_1-psm1_3) * (ph0_2-phm2_2) 
	+ 25.0*(ps0_2-psm2_2) * (phm1_1-phm1_3) - 200.0*phm1_2 + 100.0*ph0_2 + 
	100.0*phm2_2 - 200.0*phm1_2 + 100.0*phm1_1 + 100.0*phm1_3)^2 + 
	(-200.0*om0_2 + 100.0*om1_2 + 100.0*omm1_2 - 200.0*om0_2 + 100.0*om0_1 + 
	100.0*om0_3 + 1.5308084989341916e-13*ph1_2 - 1.5308084989341916e-13*phm1_2 - 
	2500.0*ph0_1 + 2500.0*ph0_3)^2 + (-200.0*ps0_2 + 100.0*ps1_2 + 
	100.0*psm1_2 - 200.0*ps0_2 + 100.0*ps0_1 + 100.0*ps0_3 + 0.25*om0_2)^2 + 
	(-25.0*(ps0_1-ps0_3) * (ph1_2-phm1_2) + 25.0*(ps1_2-psm1_2) * 
	(ph0_1-ph0_3) - 200.0*ph0_2 + 100.0*ph1_2 + 100.0*phm1_2 - 200.0*ph0_2 + 
	100.0*ph0_1 + 100.0*ph0_3)^2 + (-200.0*om1_2 + 100.0*om2_2 + 100.0*om0_2 - 
	200.0*om1_2 + 100.0*om1_1 + 100.0*om1_3 + 1.5308084989341916e-13*ph2_2 - 
	1.5308084989341916e-13*ph0_2 - 2500.0*ph1_1 + 2500.0*ph1_3)^2 + 
	(-200.0*ps1_2 + 100.0*ps2_2 + 100.0*ps0_2 - 200.0*ps1_2 + 100.0*ps1_1 + 
	100.0*ps1_3 + 0.25*om1_2)^2 + (-25.0*(ps1_1-ps1_3) * (ph2_2-ph0_2) + 
	25.0*(ps2_2-ps0_2) * (ph1_1-ph1_3) - 200.0*ph1_2 + 100.0*ph2_2 + 
	100.0*ph0_2 - 200.0*ph1_2 + 100.0*ph1_1 + 100.0*ph1_3)^2 + (-200.0*om2_2 + 
	100.0*om3_2 + 100.0*om1_2 - 200.0*om2_2 + 100.0*om2_1 + 100.0*om2_3 + 
	1.5308084989341916e-13*ph3_2 - 1.5308084989341916e-13*ph1_2 - 2500.0*ph2_1 + 
	2500.0*ph2_3)^2 + (-200.0*ps2_2 + 100.0*ps3_2 + 100.0*ps1_2 - 200.0*ps2_2 
	+ 100.0*ps2_1 + 100.0*ps2_3 + 0.25*om2_2)^2 + (-25.0*(ps2_1-ps2_3) * 
	(ph3_2-ph1_2) + 25.0*(ps3_2-ps1_2) * (ph2_1-ph2_3) - 200.0*ph2_2 + 
	100.0*ph3_2 + 100.0*ph1_2 - 200.0*ph2_2 + 100.0*ph2_1 + 100.0*ph2_3)^2 + 
	(-200.0*om3_2 + 100.0*om4_2 + 100.0*om2_2 - 200.0*om3_2 + 100.0*om3_1 + 
	100.0*om3_3 + 1.5308084989341916e-13*ph4_2 - 1.5308084989341916e-13*ph2_2 - 
	2500.0*ph3_1 + 2500.0*ph3_3)^2 + (-200.0*ps3_2 + 100.0*ps4_2 + 100.0*ps2_2 
	- 200.0*ps3_2 + 100.0*ps3_1 + 100.0*ps3_3 + 0.25*om3_2)^2 + 
	(-25.0*(ps3_1-ps3_3) * (ph4_2-ph2_2) + 25.0*(ps4_2-ps2_2) * 
	(ph3_1-ph3_3) - 200.0*ph3_2 + 100.0*ph4_2 + 100.0*ph2_2 - 200.0*ph3_2 + 
	100.0*ph3_1 + 100.0*ph3_3)^2 + (-200.0*om4_2 + 100.0*om5_2 + 100.0*om3_2 - 
	200.0*om4_2 + 100.0*om4_1 + 100.0*om4_3 + 1.5308084989341916e-13*ph5_2 - 
	1.5308084989341916e-13*ph3_2 - 2500.0*ph4_1 + 2500.0*ph4_3)^2 + 
	(-200.0*ps4_2 + 100.0*ps5_2 + 100.0*ps3_2 - 200.0*ps4_2 + 100.0*ps4_1 + 
	100.0*ps4_3 + 0.25*om4_2)^2 + (-25.0*(ps4_1-ps4_3) * (ph5_2-ph3_2) + 
	25.0*(ps5_2-ps3_2) * (ph4_1-ph4_3) - 200.0*ph4_2 + 100.0*ph5_2 + 
	100.0*ph3_2 - 200.0*ph4_2 + 100.0*ph4_1 + 100.0*ph4_3)^2 + (-200.0*om5_2 + 
	100.0*om6_2 + 100.0*om4_2 - 200.0*om5_2 + 100.0*om5_1 + 100.0*om5_3 + 
	1.5308084989341916e-13*ph6_2 - 1.5308084989341916e-13*ph4_2 - 2500.0*ph5_1 + 
	2500.0*ph5_3)^2 + (-200.0*ps5_2 + 100.0*ps6_2 + 100.0*ps4_2 - 200.0*ps5_2 
	+ 100.0*ps5_1 + 100.0*ps5_3 + 0.25*om5_2)^2 + (-25.0*(ps5_1-ps5_3) * 
	(ph6_2-ph4_2) + 25.0*(ps6_2-ps4_2) * (ph5_1-ph5_3) - 200.0*ph5_2 + 
	100.0*ph6_2 + 100.0*ph4_2 - 200.0*ph5_2 + 100.0*ph5_1 + 100.0*ph5_3)^2 + 
	(-200.0*om6_2 + 100.0*om7_2 + 100.0*om5_2 - 200.0*om6_2 + 100.0*om6_1 + 
	100.0*om6_3 + 1.5308084989341916e-13*ph7_2 - 1.5308084989341916e-13*ph5_2 - 
	2500.0*ph6_1 + 2500.0*ph6_3)^2 + (-200.0*ps6_2 + 100.0*ps7_2 + 100.0*ps5_2 
	- 200.0*ps6_2 + 100.0*ps6_1 + 100.0*ps6_3 + 0.25*om6_2)^2 + 
	(-25.0*(ps6_1-ps6_3) * (ph7_2-ph5_2) + 25.0*(ps7_2-ps5_2) * 
	(ph6_1-ph6_3) - 200.0*ph6_2 + 100.0*ph7_2 + 100.0*ph5_2 - 200.0*ph6_2 + 
	100.0*ph6_1 + 100.0*ph6_3)^2 + (-200.0*om7_2 + 100.0*om8_2 + 100.0*om6_2 - 
	200.0*om7_2 + 100.0*om7_1 + 100.0*om7_3 + 1.5308084989341916e-13*ph8_2 - 
	1.5308084989341916e-13*ph6_2 - 2500.0*ph7_1 + 2500.0*ph7_3)^2 + 
	(-200.0*ps7_2 + 100.0*ps8_2 + 100.0*ps6_2 - 200.0*ps7_2 + 100.0*ps7_1 + 
	100.0*ps7_3 + 0.25*om7_2)^2 + (-25.0*(ps7_1-ps7_3) * (ph8_2-ph6_2) + 
	25.0*(ps8_2-ps6_2) * (ph7_1-ph7_3) - 200.0*ph7_2 + 100.0*ph8_2 + 
	100.0*ph6_2 - 200.0*ph7_2 + 100.0*ph7_1 + 100.0*ph7_3)^2 + (-200.0*om8_2 + 
	100.0*om9_2 + 100.0*om7_2 - 200.0*om8_2 + 100.0*om8_1 + 100.0*om8_3 + 
	1.5308084989341916e-13*ph9_2 - 1.5308084989341916e-13*ph7_2 - 2500.0*ph8_1 + 
	2500.0*ph8_3)^2 + (-200.0*ps8_2 + 100.0*ps9_2 + 100.0*ps7_2 - 200.0*ps8_2 
	+ 100.0*ps8_1 + 100.0*ps8_3 + 0.25*om8_2)^2 + (-25.0*(ps8_1-ps8_3) * 
	(ph9_2-ph7_2) + 25.0*(ps9_2-ps7_2) * (ph8_1-ph8_3) - 200.0*ph8_2 + 
	100.0*ph9_2 + 100.0*ph7_2 - 200.0*ph8_2 + 100.0*ph8_1 + 100.0*ph8_3)^2 + 
	(-200.0*om9_2 + 100.0*om10_2 + 100.0*om8_2 - 200.0*om9_2 + 100.0*om9_1 + 
	100.0*om9_3 + 1.5308084989341916e-13*ph10_2 - 1.5308084989341916e-13*ph8_2 - 
	2500.0*ph9_1 + 2500.0*ph9_3)^2 + (-200.0*ps9_2 + 100.0*ps10_2 + 
	100.0*ps8_2 - 200.0*ps9_2 + 100.0*ps9_1 + 100.0*ps9_3 + 0.25*om9_2)^2 + 
	(-25.0*(ps9_1-ps9_3) * (ph10_2-ph8_2) + 25.0*(ps10_2-ps8_2) * 
	(ph9_1-ph9_3) - 200.0*ph9_2 + 100.0*ph10_2 + 100.0*ph8_2 - 200.0*ph9_2 + 
	100.0*ph9_1 + 100.0*ph9_3)^2 + (-200.0*omm9_1 + 100.0*omm8_1 + 
	100.0*omm10_1 - 200.0*omm9_1 + 100.0*omm9_0 + 100.0*omm9_2 + 
	1.5308084989341916e-13*phm8_1 - 1.5308084989341916e-13*phm10_1 - 
	2500.0*phm9_0 + 2500.0*phm9_2)^2 + (-200.0*psm9_1 + 100.0*psm8_1 + 
	100.0*psm10_1 - 200.0*psm9_1 + 100.0*psm9_0 + 100.0*psm9_2 + 0.25*omm9_1)^2 
	+ (-25.0*(psm9_0-psm9_2) * (phm8_1-phm10_1) + 25.0*(psm8_1-psm10_1) * 
	(phm9_0-phm9_2) - 200.0*phm9_1 + 100.0*phm8_1 + 100.0*phm10_1 - 
	200.0*phm9_1 + 100.0*phm9_0 + 100.0*phm9_2)^2 + (-200.0*omm8_1 + 
	100.0*omm7_1 + 100.0*omm9_1 - 200.0*omm8_1 + 100.0*omm8_0 + 100.0*omm8_2 + 
	1.5308084989341916e-13*phm7_1 - 1.5308084989341916e-13*phm9_1 - 2500.0*phm8_0 
	+ 2500.0*phm8_2)^2 + (-200.0*psm8_1 + 100.0*psm7_1 + 100.0*psm9_1 - 
	200.0*psm8_1 + 100.0*psm8_0 + 100.0*psm8_2 + 0.25*omm8_1)^2 + 
	(-25.0*(psm8_0-psm8_2) * (phm7_1-phm9_1) + 25.0*(psm7_1-psm9_1) * 
	(phm8_0-phm8_2) - 200.0*phm8_1 + 100.0*phm7_1 + 100.0*phm9_1 - 
	200.0*phm8_1 + 100.0*phm8_0 + 100.0*phm8_2)^2 + (-200.0*omm7_1 + 
	100.0*omm6_1 + 100.0*omm8_1 - 200.0*omm7_1 + 100.0*omm7_0 + 100.0*omm7_2 + 
	1.5308084989341916e-13*phm6_1 - 1.5308084989341916e-13*phm8_1 - 2500.0*phm7_0 
	+ 2500.0*phm7_2)^2 + (-200.0*psm7_1 + 100.0*psm6_1 + 100.0*psm8_1 - 
	200.0*psm7_1 + 100.0*psm7_0 + 100.0*psm7_2 + 0.25*omm7_1)^2 + 
	(-25.0*(psm7_0-psm7_2) * (phm6_1-phm8_1) + 25.0*(psm6_1-psm8_1) * 
	(phm7_0-phm7_2) - 200.0*phm7_1 + 100.0*phm6_1 + 100.0*phm8_1 - 
	200.0*phm7_1 + 100.0*phm7_0 + 100.0*phm7_2)^2 + (-200.0*omm6_1 + 
	100.0*omm5_1 + 100.0*omm7_1 - 200.0*omm6_1 + 100.0*omm6_0 + 100.0*omm6_2 + 
	1.5308084989341916e-13*phm5_1 - 1.5308084989341916e-13*phm7_1 - 2500.0*phm6_0 
	+ 2500.0*phm6_2)^2 + (-200.0*psm6_1 + 100.0*psm5_1 + 100.0*psm7_1 - 
	200.0*psm6_1 + 100.0*psm6_0 + 100.0*psm6_2 + 0.25*omm6_1)^2 + 
	(-25.0*(psm6_0-psm6_2) * (phm5_1-phm7_1) + 25.0*(psm5_1-psm7_1) * 
	(phm6_0-phm6_2) - 200.0*phm6_1 + 100.0*phm5_1 + 100.0*phm7_1 - 
	200.0*phm6_1 + 100.0*phm6_0 + 100.0*phm6_2)^2 + (-200.0*omm5_1 + 
	100.0*omm4_1 + 100.0*omm6_1 - 200.0*omm5_1 + 100.0*omm5_0 + 100.0*omm5_2 + 
	1.5308084989341916e-13*phm4_1 - 1.5308084989341916e-13*phm6_1 - 2500.0*phm5_0 
	+ 2500.0*phm5_2)^2 + (-200.0*psm5_1 + 100.0*psm4_1 + 100.0*psm6_1 - 
	200.0*psm5_1 + 100.0*psm5_0 + 100.0*psm5_2 + 0.25*omm5_1)^2 + 
	(-25.0*(psm5_0-psm5_2) * (phm4_1-phm6_1) + 25.0*(psm4_1-psm6_1) * 
	(phm5_0-phm5_2) - 200.0*phm5_1 + 100.0*phm4_1 + 100.0*phm6_1 - 
	200.0*phm5_1 + 100.0*phm5_0 + 100.0*phm5_2)^2 + (-200.0*omm4_1 + 
	100.0*omm3_1 + 100.0*omm5_1 - 200.0*omm4_1 + 100.0*omm4_0 + 100.0*omm4_2 + 
	1.5308084989341916e-13*phm3_1 - 1.5308084989341916e-13*phm5_1 - 2500.0*phm4_0 
	+ 2500.0*phm4_2)^2 + (-200.0*psm4_1 + 100.0*psm3_1 + 100.0*psm5_1 - 
	200.0*psm4_1 + 100.0*psm4_0 + 100.0*psm4_2 + 0.25*omm4_1)^2 + 
	(-25.0*(psm4_0-psm4_2) * (phm3_1-phm5_1) + 25.0*(psm3_1-psm5_1) * 
	(phm4_0-phm4_2) - 200.0*phm4_1 + 100.0*phm3_1 + 100.0*phm5_1 - 
	200.0*phm4_1 + 100.0*phm4_0 + 100.0*phm4_2)^2 + (-200.0*omm3_1 + 
	100.0*omm2_1 + 100.0*omm4_1 - 200.0*omm3_1 + 100.0*omm3_0 + 100.0*omm3_2 + 
	1.5308084989341916e-13*phm2_1 - 1.5308084989341916e-13*phm4_1 - 2500.0*phm3_0 
	+ 2500.0*phm3_2)^2 + (-200.0*psm3_1 + 100.0*psm2_1 + 100.0*psm4_1 - 
	200.0*psm3_1 + 100.0*psm3_0 + 100.0*psm3_2 + 0.25*omm3_1)^2 + 
	(-25.0*(psm3_0-psm3_2) * (phm2_1-phm4_1) + 25.0*(psm2_1-psm4_1) * 
	(phm3_0-phm3_2) - 200.0*phm3_1 + 100.0*phm2_1 + 100.0*phm4_1 - 
	200.0*phm3_1 + 100.0*phm3_0 + 100.0*phm3_2)^2 + (-200.0*omm2_1 + 
	100.0*omm1_1 + 100.0*omm3_1 - 200.0*omm2_1 + 100.0*omm2_0 + 100.0*omm2_2 + 
	1.5308084989341916e-13*phm1_1 - 1.5308084989341916e-13*phm3_1 - 2500.0*phm2_0 
	+ 2500.0*phm2_2)^2 + (-200.0*psm2_1 + 100.0*psm1_1 + 100.0*psm3_1 - 
	200.0*psm2_1 + 100.0*psm2_0 + 100.0*psm2_2 + 0.25*omm2_1)^2 + 
	(-25.0*(psm2_0-psm2_2) * (phm1_1-phm3_1) + 25.0*(psm1_1-psm3_1) * 
	(phm2_0-phm2_2) - 200.0*phm2_1 + 100.0*phm1_1 + 100.0*phm3_1 - 
	200.0*phm2_1 + 100.0*phm2_0 + 100.0*phm2_2)^2 + (-200.0*omm1_1 + 
	100.0*om0_1 + 100.0*omm2_1 - 200.0*omm1_1 + 100.0*omm1_0 + 100.0*omm1_2 + 
	1.5308084989341916e-13*ph0_1 - 1.5308084989341916e-13*phm2_1 - 2500.0*phm1_0 
	+ 2500.0*phm1_2)^2 + (-200.0*psm1_1 + 100.0*ps0_1 + 100.0*psm2_1 - 
	200.0*psm1_1 + 100.0*psm1_0 + 100.0*psm1_2 + 0.25*omm1_1)^2 + 
	(-25.0*(psm1_0-psm1_2) * (ph0_1-phm2_1) + 25.0*(ps0_1-psm2_1) * 
	(phm1_0-phm1_2) - 200.0*phm1_1 + 100.0*ph0_1 + 100.0*phm2_1 - 200.0*phm1_1 
	+ 100.0*phm1_0 + 100.0*phm1_2)^2 + (-200.0*om0_1 + 100.0*om1_1 + 
	100.0*omm1_1 - 200.0*om0_1 + 100.0*om0_0 + 100.0*om0_2 + 
	1.5308084989341916e-13*ph1_1 - 1.5308084989341916e-13*phm1_1 - 2500.0*ph0_0 + 
	2500.0*ph0_2)^2 + (-200.0*ps0_1 + 100.0*ps1_1 + 100.0*psm1_1 - 200.0*ps0_1 
	+ 100.0*ps0_0 + 100.0*ps0_2 + 0.25*om0_1)^2 + (-25.0*(ps0_0-ps0_2) * 
	(ph1_1-phm1_1) + 25.0*(ps1_1-psm1_1) * (ph0_0-ph0_2) - 200.0*ph0_1 + 
	100.0*ph1_1 + 100.0*phm1_1 - 200.0*ph0_1 + 100.0*ph0_0 + 100.0*ph0_2)^2 + 
	(-200.0*om1_1 + 100.0*om2_1 + 100.0*om0_1 - 200.0*om1_1 + 100.0*om1_0 + 
	100.0*om1_2 + 1.5308084989341916e-13*ph2_1 - 1.5308084989341916e-13*ph0_1 - 
	2500.0*ph1_0 + 2500.0*ph1_2)^2 + (-200.0*ps1_1 + 100.0*ps2_1 + 100.0*ps0_1 
	- 200.0*ps1_1 + 100.0*ps1_0 + 100.0*ps1_2 + 0.25*om1_1)^2 + 
	(-25.0*(ps1_0-ps1_2) * (ph2_1-ph0_1) + 25.0*(ps2_1-ps0_1) * (ph1_0-ph1_2) 
	- 200.0*ph1_1 + 100.0*ph2_1 + 100.0*ph0_1 - 200.0*ph1_1 + 100.0*ph1_0 + 
	100.0*ph1_2)^2 + (-200.0*om2_1 + 100.0*om3_1 + 100.0*om1_1 - 200.0*om2_1 + 
	100.0*om2_0 + 100.0*om2_2 + 1.5308084989341916e-13*ph3_1 - 
	1.5308084989341916e-13*ph1_1 - 2500.0*ph2_0 + 2500.0*ph2_2)^2 + 
	(-200.0*ps2_1 + 100.0*ps3_1 + 100.0*ps1_1 - 200.0*ps2_1 + 100.0*ps2_0 + 
	100.0*ps2_2 + 0.25*om2_1)^2 + (-25.0*(ps2_0-ps2_2) * (ph3_1-ph1_1) + 
	25.0*(ps3_1-ps1_1) * (ph2_0-ph2_2) - 200.0*ph2_1 + 100.0*ph3_1 + 
	100.0*ph1_1 - 200.0*ph2_1 + 100.0*ph2_0 + 100.0*ph2_2)^2 + (-200.0*om3_1 + 
	100.0*om4_1 + 100.0*om2_1 - 200.0*om3_1 + 100.0*om3_0 + 100.0*om3_2 + 
	1.5308084989341916e-13*ph4_1 - 1.5308084989341916e-13*ph2_1 - 2500.0*ph3_0 + 
	2500.0*ph3_2)^2 + (-200.0*ps3_1 + 100.0*ps4_1 + 100.0*ps2_1 - 200.0*ps3_1 
	+ 100.0*ps3_0 + 100.0*ps3_2 + 0.25*om3_1)^2 + (-25.0*(ps3_0-ps3_2) * 
	(ph4_1-ph2_1) + 25.0*(ps4_1-ps2_1) * (ph3_0-ph3_2) - 200.0*ph3_1 + 
	100.0*ph4_1 + 100.0*ph2_1 - 200.0*ph3_1 + 100.0*ph3_0 + 100.0*ph3_2)^2 + 
	(-200.0*om4_1 + 100.0*om5_1 + 100.0*om3_1 - 200.0*om4_1 + 100.0*om4_0 + 
	100.0*om4_2 + 1.5308084989341916e-13*ph5_1 - 1.5308084989341916e-13*ph3_1 - 
	2500.0*ph4_0 + 2500.0*ph4_2)^2 + (-200.0*ps4_1 + 100.0*ps5_1 + 100.0*ps3_1 
	- 200.0*ps4_1 + 100.0*ps4_0 + 100.0*ps4_2 + 0.25*om4_1)^2 + 
	(-25.0*(ps4_0-ps4_2) * (ph5_1-ph3_1) + 25.0*(ps5_1-ps3_1) * (ph4_0-ph4_2) 
	- 200.0*ph4_1 + 100.0*ph5_1 + 100.0*ph3_1 - 200.0*ph4_1 + 100.0*ph4_0 + 
	100.0*ph4_2)^2 + (-200.0*om5_1 + 100.0*om6_1 + 100.0*om4_1 - 200.0*om5_1 + 
	100.0*om5_0 + 100.0*om5_2 + 1.5308084989341916e-13*ph6_1 - 
	1.5308084989341916e-13*ph4_1 - 2500.0*ph5_0 + 2500.0*ph5_2)^2 + 
	(-200.0*ps5_1 + 100.0*ps6_1 + 100.0*ps4_1 - 200.0*ps5_1 + 100.0*ps5_0 + 
	100.0*ps5_2 + 0.25*om5_1)^2 + (-25.0*(ps5_0-ps5_2) * (ph6_1-ph4_1) + 
	25.0*(ps6_1-ps4_1) * (ph5_0-ph5_2) - 200.0*ph5_1 + 100.0*ph6_1 + 
	100.0*ph4_1 - 200.0*ph5_1 + 100.0*ph5_0 + 100.0*ph5_2)^2 + (-200.0*om6_1 + 
	100.0*om7_1 + 100.0*om5_1 - 200.0*om6_1 + 100.0*om6_0 + 100.0*om6_2 + 
	1.5308084989341916e-13*ph7_1 - 1.5308084989341916e-13*ph5_1 - 2500.0*ph6_0 + 
	2500.0*ph6_2)^2 + (-200.0*ps6_1 + 100.0*ps7_1 + 100.0*ps5_1 - 200.0*ps6_1 
	+ 100.0*ps6_0 + 100.0*ps6_2 + 0.25*om6_1)^2 + (-25.0*(ps6_0-ps6_2) * 
	(ph7_1-ph5_1) + 25.0*(ps7_1-ps5_1) * (ph6_0-ph6_2) - 200.0*ph6_1 + 
	100.0*ph7_1 + 100.0*ph5_1 - 200.0*ph6_1 + 100.0*ph6_0 + 100.0*ph6_2)^2 + 
	(-200.0*om7_1 + 100.0*om8_1 + 100.0*om6_1 - 200.0*om7_1 + 100.0*om7_0 + 
	100.0*om7_2 + 1.5308084989341916e-13*ph8_1 - 1.5308084989341916e-13*ph6_1 - 
	2500.0*ph7_0 + 2500.0*ph7_2)^2 + (-200.0*ps7_1 + 100.0*ps8_1 + 100.0*ps6_1 
	- 200.0*ps7_1 + 100.0*ps7_0 + 100.0*ps7_2 + 0.25*om7_1)^2 + 
	(-25.0*(ps7_0-ps7_2) * (ph8_1-ph6_1) + 25.0*(ps8_1-ps6_1) * (ph7_0-ph7_2) 
	- 200.0*ph7_1 + 100.0*ph8_1 + 100.0*ph6_1 - 200.0*ph7_1 + 100.0*ph7_0 + 
	100.0*ph7_2)^2 + (-200.0*om8_1 + 100.0*om9_1 + 100.0*om7_1 - 200.0*om8_1 + 
	100.0*om8_0 + 100.0*om8_2 + 1.5308084989341916e-13*ph9_1 - 
	1.5308084989341916e-13*ph7_1 - 2500.0*ph8_0 + 2500.0*ph8_2)^2 + 
	(-200.0*ps8_1 + 100.0*ps9_1 + 100.0*ps7_1 - 200.0*ps8_1 + 100.0*ps8_0 + 
	100.0*ps8_2 + 0.25*om8_1)^2 + (-25.0*(ps8_0-ps8_2) * (ph9_1-ph7_1) + 
	25.0*(ps9_1-ps7_1) * (ph8_0-ph8_2) - 200.0*ph8_1 + 100.0*ph9_1 + 
	100.0*ph7_1 - 200.0*ph8_1 + 100.0*ph8_0 + 100.0*ph8_2)^2 + (-200.0*om9_1 + 
	100.0*om10_1 + 100.0*om8_1 - 200.0*om9_1 + 100.0*om9_0 + 100.0*om9_2 + 
	1.5308084989341916e-13*ph10_1 - 1.5308084989341916e-13*ph8_1 - 2500.0*ph9_0 + 
	2500.0*ph9_2)^2 + (-200.0*ps9_1 + 100.0*ps10_1 + 100.0*ps8_1 - 200.0*ps9_1 
	+ 100.0*ps9_0 + 100.0*ps9_2 + 0.25*om9_1)^2 + (-25.0*(ps9_0-ps9_2) * 
	(ph10_1-ph8_1) + 25.0*(ps10_1-ps8_1) * (ph9_0-ph9_2) - 200.0*ph9_1 + 
	100.0*ph10_1 + 100.0*ph8_1 - 200.0*ph9_1 + 100.0*ph9_0 + 100.0*ph9_2)^2 + 
	(-200.0*omm9_0 + 100.0*omm8_0 + 100.0*omm10_0 - 200.0*omm9_0 + 100.0*omm9_1 + 
	100.0*omm9_1 + 1.5308084989341916e-13*phm8_0 - 1.5308084989341916e-13*phm10_0 
	- 2500.0*phm9_1 + 2500.0*phm9_1)^2 + (-200.0*psm9_0 + 100.0*psm8_0 + 
	100.0*psm10_0 - 200.0*psm9_0 + 100.0*psm9_1 + 100.0*psm9_1 + 0.25*omm9_0)^2 + 
	(-25.0*(psm9_1-psm9_1) * (phm8_0-phm10_0) + 25.0*(psm8_0-psm10_0) * 
	(phm9_1-phm9_1) - 200.0*phm9_0 + 100.0*phm8_0 + 100.0*phm10_0 - 200.0*phm9_0 + 
	100.0*phm9_1 + 100.0*phm9_1)^2 + (-200.0*omm8_0 + 100.0*omm7_0 + 100.0*omm9_0 
	- 200.0*omm8_0 + 100.0*omm8_1 + 100.0*omm8_1 + 1.5308084989341916e-13*phm7_0 - 
	1.5308084989341916e-13*phm9_0 - 2500.0*phm8_1 + 2500.0*phm8_1)^2 + 
	(-200.0*psm8_0 + 100.0*psm7_0 + 100.0*psm9_0 - 200.0*psm8_0 + 100.0*psm8_1 + 
	100.0*psm8_1 + 0.25*omm8_0)^2 + (-25.0*(psm8_1-psm8_1) * (phm7_0-phm9_0) + 
	25.0*(psm7_0-psm9_0) * (phm8_1-phm8_1) - 200.0*phm8_0 + 100.0*phm7_0 + 
	100.0*phm9_0 - 200.0*phm8_0 + 100.0*phm8_1 + 100.0*phm8_1)^2 + (-200.0*omm7_0 
	+ 100.0*omm6_0 + 100.0*omm8_0 - 200.0*omm7_0 + 100.0*omm7_1 + 100.0*omm7_1 + 
	1.5308084989341916e-13*phm6_0 - 1.5308084989341916e-13*phm8_0 - 2500.0*phm7_1 + 
	2500.0*phm7_1)^2 + (-200.0*psm7_0 + 100.0*psm6_0 + 100.0*psm8_0 - 200.0*psm7_0 
	+ 100.0*psm7_1 + 100.0*psm7_1 + 0.25*omm7_0)^2 + (-25.0*(psm7_1-psm7_1) * 
	(phm6_0-phm8_0) + 25.0*(psm6_0-psm8_0) * (phm7_1-phm7_1) - 200.0*phm7_0 + 
	100.0*phm6_0 + 100.0*phm8_0 - 200.0*phm7_0 + 100.0*phm7_1 + 100.0*phm7_1)^2 + 
	(-200.0*omm6_0 + 100.0*omm5_0 + 100.0*omm7_0 - 200.0*omm6_0 + 100.0*omm6_1 + 
	100.0*omm6_1 + 1.5308084989341916e-13*phm5_0 - 1.5308084989341916e-13*phm7_0 - 
	2500.0*phm6_1 + 2500.0*phm6_1)^2 + (-200.0*psm6_0 + 100.0*psm5_0 + 
	100.0*psm7_0 - 200.0*psm6_0 + 100.0*psm6_1 + 100.0*psm6_1 + 0.25*omm6_0)^2 + 
	(-25.0*(psm6_1-psm6_1) * (phm5_0-phm7_0) + 25.0*(psm5_0-psm7_0) * 
	(phm6_1-phm6_1) - 200.0*phm6_0 + 100.0*phm5_0 + 100.0*phm7_0 - 200.0*phm6_0 + 
	100.0*phm6_1 + 100.0*phm6_1)^2 + (-200.0*omm5_0 + 100.0*omm4_0 + 100.0*omm6_0 
	- 200.0*omm5_0 + 100.0*omm5_1 + 100.0*omm5_1 + 1.5308084989341916e-13*phm4_0 - 
	1.5308084989341916e-13*phm6_0 - 2500.0*phm5_1 + 2500.0*phm5_1)^2 + 
	(-200.0*psm5_0 + 100.0*psm4_0 + 100.0*psm6_0 - 200.0*psm5_0 + 100.0*psm5_1 + 
	100.0*psm5_1 + 0.25*omm5_0)^2 + (-25.0*(psm5_1-psm5_1) * (phm4_0-phm6_0) + 
	25.0*(psm4_0-psm6_0) * (phm5_1-phm5_1) - 200.0*phm5_0 + 100.0*phm4_0 + 
	100.0*phm6_0 - 200.0*phm5_0 + 100.0*phm5_1 + 100.0*phm5_1)^2 + (-200.0*omm4_0 
	+ 100.0*omm3_0 + 100.0*omm5_0 - 200.0*omm4_0 + 100.0*omm4_1 + 100.0*omm4_1 + 
	1.5308084989341916e-13*phm3_0 - 1.5308084989341916e-13*phm5_0 - 2500.0*phm4_1 + 
	2500.0*phm4_1)^2 + (-200.0*psm4_0 + 100.0*psm3_0 + 100.0*psm5_0 - 200.0*psm4_0 
	+ 100.0*psm4_1 + 100.0*psm4_1 + 0.25*omm4_0)^2 + (-25.0*(psm4_1-psm4_1) * 
	(phm3_0-phm5_0) + 25.0*(psm3_0-psm5_0) * (phm4_1-phm4_1) - 200.0*phm4_0 + 
	100.0*phm3_0 + 100.0*phm5_0 - 200.0*phm4_0 + 100.0*phm4_1 + 100.0*phm4_1)^2 + 
	(-200.0*omm3_0 + 100.0*omm2_0 + 100.0*omm4_0 - 200.0*omm3_0 + 100.0*omm3_1 + 
	100.0*omm3_1 + 1.5308084989341916e-13*phm2_0 - 1.5308084989341916e-13*phm4_0 - 
	2500.0*phm3_1 + 2500.0*phm3_1)^2 + (-200.0*psm3_0 + 100.0*psm2_0 + 
	100.0*psm4_0 - 200.0*psm3_0 + 100.0*psm3_1 + 100.0*psm3_1 + 0.25*omm3_0)^2 + 
	(-25.0*(psm3_1-psm3_1) * (phm2_0-phm4_0) + 25.0*(psm2_0-psm4_0) * 
	(phm3_1-phm3_1) - 200.0*phm3_0 + 100.0*phm2_0 + 100.0*phm4_0 - 200.0*phm3_0 + 
	100.0*phm3_1 + 100.0*phm3_1)^2 + (-200.0*omm2_0 + 100.0*omm1_0 + 100.0*omm3_0 
	- 200.0*omm2_0 + 100.0*omm2_1 + 100.0*omm2_1 + 1.5308084989341916e-13*phm1_0 - 
	1.5308084989341916e-13*phm3_0 - 2500.0*phm2_1 + 2500.0*phm2_1)^2 + 
	(-200.0*psm2_0 + 100.0*psm1_0 + 100.0*psm3_0 - 200.0*psm2_0 + 100.0*psm2_1 + 
	100.0*psm2_1 + 0.25*omm2_0)^2 + (-25.0*(psm2_1-psm2_1) * (phm1_0-phm3_0) + 
	25.0*(psm1_0-psm3_0) * (phm2_1-phm2_1) - 200.0*phm2_0 + 100.0*phm1_0 + 
	100.0*phm3_0 - 200.0*phm2_0 + 100.0*phm2_1 + 100.0*phm2_1)^2 + (-200.0*omm1_0 
	+ 100.0*om0_0 + 100.0*omm2_0 - 200.0*omm1_0 + 100.0*omm1_1 + 100.0*omm1_1 + 
	1.5308084989341916e-13*ph0_0 - 1.5308084989341916e-13*phm2_0 - 2500.0*phm1_1 + 
	2500.0*phm1_1)^2 + (-200.0*psm1_0 + 100.0*ps0_0 + 100.0*psm2_0 - 200.0*psm1_0 
	+ 100.0*psm1_1 + 100.0*psm1_1 + 0.25*omm1_0)^2 + (-25.0*(psm1_1-psm1_1) * 
	(ph0_0-phm2_0) + 25.0*(ps0_0-psm2_0) * (phm1_1-phm1_1) - 200.0*phm1_0 + 
	100.0*ph0_0 + 100.0*phm2_0 - 200.0*phm1_0 + 100.0*phm1_1 + 100.0*phm1_1)^2 + 
	(-200.0*om0_0 + 100.0*om1_0 + 100.0*omm1_0 - 200.0*om0_0 + 100.0*om0_1 + 
	100.0*om0_1 + 1.5308084989341916e-13*ph1_0 - 1.5308084989341916e-13*phm1_0 - 
	2500.0*ph0_1 + 2500.0*ph0_1)^2 + (-200.0*ps0_0 + 100.0*ps1_0 + 100.0*psm1_0 - 
	200.0*ps0_0 + 100.0*ps0_1 + 100.0*ps0_1 + 0.25*om0_0)^2 + 
	(-25.0*(ps0_1-ps0_1) * (ph1_0-phm1_0) + 25.0*(ps1_0-psm1_0) * (ph0_1-ph0_1) - 
	200.0*ph0_0 + 100.0*ph1_0 + 100.0*phm1_0 - 200.0*ph0_0 + 100.0*ph0_1 + 
	100.0*ph0_1)^2 + (-200.0*om1_0 + 100.0*om2_0 + 100.0*om0_0 - 200.0*om1_0 + 
	100.0*om1_1 + 100.0*om1_1 + 1.5308084989341916e-13*ph2_0 - 
	1.5308084989341916e-13*ph0_0 - 2500.0*ph1_1 + 2500.0*ph1_1)^2 + (-200.0*ps1_0 
	+ 100.0*ps2_0 + 100.0*ps0_0 - 200.0*ps1_0 + 100.0*ps1_1 + 100.0*ps1_1 + 
	0.25*om1_0)^2 + (-25.0*(ps1_1-ps1_1) * (ph2_0-ph0_0) + 25.0*(ps2_0-ps0_0) * 
	(ph1_1-ph1_1) - 200.0*ph1_0 + 100.0*ph2_0 + 100.0*ph0_0 - 200.0*ph1_0 + 
	100.0*ph1_1 + 100.0*ph1_1)^2 + (-200.0*om2_0 + 100.0*om3_0 + 100.0*om1_0 - 
	200.0*om2_0 + 100.0*om2_1 + 100.0*om2_1 + 1.5308084989341916e-13*ph3_0 - 
	1.5308084989341916e-13*ph1_0 - 2500.0*ph2_1 + 2500.0*ph2_1)^2 + (-200.0*ps2_0 
	+ 100.0*ps3_0 + 100.0*ps1_0 - 200.0*ps2_0 + 100.0*ps2_1 + 100.0*ps2_1 + 
	0.25*om2_0)^2 + (-25.0*(ps2_1-ps2_1) * (ph3_0-ph1_0) + 25.0*(ps3_0-ps1_0) * 
	(ph2_1-ph2_1) - 200.0*ph2_0 + 100.0*ph3_0 + 100.0*ph1_0 - 200.0*ph2_0 + 
	100.0*ph2_1 + 100.0*ph2_1)^2 + (-200.0*om3_0 + 100.0*om4_0 + 100.0*om2_0 - 
	200.0*om3_0 + 100.0*om3_1 + 100.0*om3_1 + 1.5308084989341916e-13*ph4_0 - 
	1.5308084989341916e-13*ph2_0 - 2500.0*ph3_1 + 2500.0*ph3_1)^2 + (-200.0*ps3_0 
	+ 100.0*ps4_0 + 100.0*ps2_0 - 200.0*ps3_0 + 100.0*ps3_1 + 100.0*ps3_1 + 
	0.25*om3_0)^2 + (-25.0*(ps3_1-ps3_1) * (ph4_0-ph2_0) + 25.0*(ps4_0-ps2_0) * 
	(ph3_1-ph3_1) - 200.0*ph3_0 + 100.0*ph4_0 + 100.0*ph2_0 - 200.0*ph3_0 + 
	100.0*ph3_1 + 100.0*ph3_1)^2 + (-200.0*om4_0 + 100.0*om5_0 + 100.0*om3_0 - 
	200.0*om4_0 + 100.0*om4_1 + 100.0*om4_1 + 1.5308084989341916e-13*ph5_0 - 
	1.5308084989341916e-13*ph3_0 - 2500.0*ph4_1 + 2500.0*ph4_1)^2 + (-200.0*ps4_0 
	+ 100.0*ps5_0 + 100.0*ps3_0 - 200.0*ps4_0 + 100.0*ps4_1 + 100.0*ps4_1 + 
	0.25*om4_0)^2 + (-25.0*(ps4_1-ps4_1) * (ph5_0-ph3_0) + 25.0*(ps5_0-ps3_0) * 
	(ph4_1-ph4_1) - 200.0*ph4_0 + 100.0*ph5_0 + 100.0*ph3_0 - 200.0*ph4_0 + 
	100.0*ph4_1 + 100.0*ph4_1)^2 + (-200.0*om5_0 + 100.0*om6_0 + 100.0*om4_0 - 
	200.0*om5_0 + 100.0*om5_1 + 100.0*om5_1 + 1.5308084989341916e-13*ph6_0 - 
	1.5308084989341916e-13*ph4_0 - 2500.0*ph5_1 + 2500.0*ph5_1)^2 + (-200.0*ps5_0 
	+ 100.0*ps6_0 + 100.0*ps4_0 - 200.0*ps5_0 + 100.0*ps5_1 + 100.0*ps5_1 + 
	0.25*om5_0)^2 + (-25.0*(ps5_1-ps5_1) * (ph6_0-ph4_0) + 25.0*(ps6_0-ps4_0) * 
	(ph5_1-ph5_1) - 200.0*ph5_0 + 100.0*ph6_0 + 100.0*ph4_0 - 200.0*ph5_0 + 
	100.0*ph5_1 + 100.0*ph5_1)^2 + (-200.0*om6_0 + 100.0*om7_0 + 100.0*om5_0 - 
	200.0*om6_0 + 100.0*om6_1 + 100.0*om6_1 + 1.5308084989341916e-13*ph7_0 - 
	1.5308084989341916e-13*ph5_0 - 2500.0*ph6_1 + 2500.0*ph6_1)^2 + (-200.0*ps6_0 
	+ 100.0*ps7_0 + 100.0*ps5_0 - 200.0*ps6_0 + 100.0*ps6_1 + 100.0*ps6_1 + 
	0.25*om6_0)^2 + (-25.0*(ps6_1-ps6_1) * (ph7_0-ph5_0) + 25.0*(ps7_0-ps5_0) * 
	(ph6_1-ph6_1) - 200.0*ph6_0 + 100.0*ph7_0 + 100.0*ph5_0 - 200.0*ph6_0 + 
	100.0*ph6_1 + 100.0*ph6_1)^2 + (-200.0*om7_0 + 100.0*om8_0 + 100.0*om6_0 - 
	200.0*om7_0 + 100.0*om7_1 + 100.0*om7_1 + 1.5308084989341916e-13*ph8_0 - 
	1.5308084989341916e-13*ph6_0 - 2500.0*ph7_1 + 2500.0*ph7_1)^2 + (-200.0*ps7_0 
	+ 100.0*ps8_0 + 100.0*ps6_0 - 200.0*ps7_0 + 100.0*ps7_1 + 100.0*ps7_1 + 
	0.25*om7_0)^2 + (-25.0*(ps7_1-ps7_1) * (ph8_0-ph6_0) + 25.0*(ps8_0-ps6_0) * 
	(ph7_1-ph7_1) - 200.0*ph7_0 + 100.0*ph8_0 + 100.0*ph6_0 - 200.0*ph7_0 + 
	100.0*ph7_1 + 100.0*ph7_1)^2 + (-200.0*om8_0 + 100.0*om9_0 + 100.0*om7_0 - 
	200.0*om8_0 + 100.0*om8_1 + 100.0*om8_1 + 1.5308084989341916e-13*ph9_0 - 
	1.5308084989341916e-13*ph7_0 - 2500.0*ph8_1 + 2500.0*ph8_1)^2 + (-200.0*ps8_0 
	+ 100.0*ps9_0 + 100.0*ps7_0 - 200.0*ps8_0 + 100.0*ps8_1 + 100.0*ps8_1 + 
	0.25*om8_0)^2 + (-25.0*(ps8_1-ps8_1) * (ph9_0-ph7_0) + 25.0*(ps9_0-ps7_0) * 
	(ph8_1-ph8_1) - 200.0*ph8_0 + 100.0*ph9_0 + 100.0*ph7_0 - 200.0*ph8_0 + 
	100.0*ph8_1 + 100.0*ph8_1)^2 + (-200.0*om9_0 + 100.0*om10_0 + 100.0*om8_0 - 
	200.0*om9_0 + 100.0*om9_1 + 100.0*om9_1 + 1.5308084989341916e-13*ph10_0 - 
	1.5308084989341916e-13*ph8_0 - 2500.0*ph9_1 + 2500.0*ph9_1)^2 + (-200.0*ps9_0 
	+ 100.0*ps10_0 + 100.0*ps8_0 - 200.0*ps9_0 + 100.0*ps9_1 + 100.0*ps9_1 + 
	0.25*om9_0)^2 + (-25.0*(ps9_1-ps9_1) * (ph10_0-ph8_0) + 25.0*(ps10_0-ps8_0) * 
	(ph9_1-ph9_1) - 200.0*ph9_0 + 100.0*ph10_0 + 100.0*ph8_0 - 200.0*ph9_0 + 
	100.0*ph9_1 + 100.0*ph9_1)^2 + (-200.0*omm9_1 + 100.0*omm8_1 + 100.0*omm10_1 - 
	200.0*omm9_1 + 100.0*omm9_2 + 100.0*omm9_0 + 1.5308084989341916e-13*phm8_1 - 
	1.5308084989341916e-13*phm10_1 - 2500.0*phm9_2 + 2500.0*phm9_0)^2 + 
	(-200.0*psm9_1 + 100.0*psm8_1 + 100.0*psm10_1 - 200.0*psm9_1 + 100.0*psm9_2 + 
	100.0*psm9_0 + 0.25*omm9_1)^2 + (-25.0*(psm9_2-psm9_0) * (phm8_1-phm10_1) + 
	25.0*(psm8_1-psm10_1) * (phm9_2-phm9_0) - 200.0*phm9_1 + 100.0*phm8_1 + 
	100.0*phm10_1 - 200.0*phm9_1 + 100.0*phm9_2 + 100.0*phm9_0)^2 + (-200.0*omm8_1 
	+ 100.0*omm7_1 + 100.0*omm9_1 - 200.0*omm8_1 + 100.0*omm8_2 + 100.0*omm8_0 + 
	1.5308084989341916e-13*phm7_1 - 1.5308084989341916e-13*phm9_1 - 2500.0*phm8_2 + 
	2500.0*phm8_0)^2 + (-200.0*psm8_1 + 100.0*psm7_1 + 100.0*psm9_1 - 200.0*psm8_1 
	+ 100.0*psm8_2 + 100.0*psm8_0 + 0.25*omm8_1)^2 + (-25.0*(psm8_2-psm8_0) * 
	(phm7_1-phm9_1) + 25.0*(psm7_1-psm9_1) * (phm8_2-phm8_0) - 200.0*phm8_1 + 
	100.0*phm7_1 + 100.0*phm9_1 - 200.0*phm8_1 + 100.0*phm8_2 + 100.0*phm8_0)^2 + 
	(-200.0*omm7_1 + 100.0*omm6_1 + 100.0*omm8_1 - 200.0*omm7_1 + 100.0*omm7_2 + 
	100.0*omm7_0 + 1.5308084989341916e-13*phm6_1 - 1.5308084989341916e-13*phm8_1 - 
	2500.0*phm7_2 + 2500.0*phm7_0)^2 + (-200.0*psm7_1 + 100.0*psm6_1 + 100.0*psm8_1 
	- 200.0*psm7_1 + 100.0*psm7_2 + 100.0*psm7_0 + 0.25*omm7_1)^2 + 
	(-25.0*(psm7_2-psm7_0) * (phm6_1-phm8_1) + 25.0*(psm6_1-psm8_1) * 
	(phm7_2-phm7_0) - 200.0*phm7_1 + 100.0*phm6_1 + 100.0*phm8_1 - 200.0*phm7_1 + 
	100.0*phm7_2 + 100.0*phm7_0)^2 + (-200.0*omm6_1 + 100.0*omm5_1 + 100.0*omm7_1 - 
	200.0*omm6_1 + 100.0*omm6_2 + 100.0*omm6_0 + 1.5308084989341916e-13*phm5_1 - 
	1.5308084989341916e-13*phm7_1 - 2500.0*phm6_2 + 2500.0*phm6_0)^2 + 
	(-200.0*psm6_1 + 100.0*psm5_1 + 100.0*psm7_1 - 200.0*psm6_1 + 100.0*psm6_2 + 
	100.0*psm6_0 + 0.25*omm6_1)^2 + (-25.0*(psm6_2-psm6_0) * (phm5_1-phm7_1) + 
	25.0*(psm5_1-psm7_1) * (phm6_2-phm6_0) - 200.0*phm6_1 + 100.0*phm5_1 + 
	100.0*phm7_1 - 200.0*phm6_1 + 100.0*phm6_2 + 100.0*phm6_0)^2 + (-200.0*omm5_1 + 
	100.0*omm4_1 + 100.0*omm6_1 - 200.0*omm5_1 + 100.0*omm5_2 + 100.0*omm5_0 + 
	1.5308084989341916e-13*phm4_1 - 1.5308084989341916e-13*phm6_1 - 2500.0*phm5_2 + 
	2500.0*phm5_0)^2 + (-200.0*psm5_1 + 100.0*psm4_1 + 100.0*psm6_1 - 200.0*psm5_1 
	+ 100.0*psm5_2 + 100.0*psm5_0 + 0.25*omm5_1)^2 + (-25.0*(psm5_2-psm5_0) * 
	(phm4_1-phm6_1) + 25.0*(psm4_1-psm6_1) * (phm5_2-phm5_0) - 200.0*phm5_1 + 
	100.0*phm4_1 + 100.0*phm6_1 - 200.0*phm5_1 + 100.0*phm5_2 + 100.0*phm5_0)^2 + 
	(-200.0*omm4_1 + 100.0*omm3_1 + 100.0*omm5_1 - 200.0*omm4_1 + 100.0*omm4_2 + 
	100.0*omm4_0 + 1.5308084989341916e-13*phm3_1 - 1.5308084989341916e-13*phm5_1 - 
	2500.0*phm4_2 + 2500.0*phm4_0)^2 + (-200.0*psm4_1 + 100.0*psm3_1 + 100.0*psm5_1 
	- 200.0*psm4_1 + 100.0*psm4_2 + 100.0*psm4_0 + 0.25*omm4_1)^2 + 
	(-25.0*(psm4_2-psm4_0) * (phm3_1-phm5_1) + 25.0*(psm3_1-psm5_1) * 
	(phm4_2-phm4_0) - 200.0*phm4_1 + 100.0*phm3_1 + 100.0*phm5_1 - 200.0*phm4_1 + 
	100.0*phm4_2 + 100.0*phm4_0)^2 + (-200.0*omm3_1 + 100.0*omm2_1 + 100.0*omm4_1 - 
	200.0*omm3_1 + 100.0*omm3_2 + 100.0*omm3_0 + 1.5308084989341916e-13*phm2_1 - 
	1.5308084989341916e-13*phm4_1 - 2500.0*phm3_2 + 2500.0*phm3_0)^2 + 
	(-200.0*psm3_1 + 100.0*psm2_1 + 100.0*psm4_1 - 200.0*psm3_1 + 100.0*psm3_2 + 
	100.0*psm3_0 + 0.25*omm3_1)^2 + (-25.0*(psm3_2-psm3_0) * (phm2_1-phm4_1) + 
	25.0*(psm2_1-psm4_1) * (phm3_2-phm3_0) - 200.0*phm3_1 + 100.0*phm2_1 + 
	100.0*phm4_1 - 200.0*phm3_1 + 100.0*phm3_2 + 100.0*phm3_0)^2 + (-200.0*omm2_1 + 
	100.0*omm1_1 + 100.0*omm3_1 - 200.0*omm2_1 + 100.0*omm2_2 + 100.0*omm2_0 + 
	1.5308084989341916e-13*phm1_1 - 1.5308084989341916e-13*phm3_1 - 2500.0*phm2_2 + 
	2500.0*phm2_0)^2 + (-200.0*psm2_1 + 100.0*psm1_1 + 100.0*psm3_1 - 200.0*psm2_1 
	+ 100.0*psm2_2 + 100.0*psm2_0 + 0.25*omm2_1)^2 + (-25.0*(psm2_2-psm2_0) * 
	(phm1_1-phm3_1) + 25.0*(psm1_1-psm3_1) * (phm2_2-phm2_0) - 200.0*phm2_1 + 
	100.0*phm1_1 + 100.0*phm3_1 - 200.0*phm2_1 + 100.0*phm2_2 + 100.0*phm2_0)^2 + 
	(-200.0*omm1_1 + 100.0*om0_1 + 100.0*omm2_1 - 200.0*omm1_1 + 100.0*omm1_2 + 
	100.0*omm1_0 + 1.5308084989341916e-13*ph0_1 - 1.5308084989341916e-13*phm2_1 - 
	2500.0*phm1_2 + 2500.0*phm1_0)^2 + (-200.0*psm1_1 + 100.0*ps0_1 + 100.0*psm2_1 
	- 200.0*psm1_1 + 100.0*psm1_2 + 100.0*psm1_0 + 0.25*omm1_1)^2 + 
	(-25.0*(psm1_2-psm1_0) * (ph0_1-phm2_1) + 25.0*(ps0_1-psm2_1) * (phm1_2-phm1_0) 
	- 200.0*phm1_1 + 100.0*ph0_1 + 100.0*phm2_1 - 200.0*phm1_1 + 100.0*phm1_2 + 
	100.0*phm1_0)^2 + (-200.0*om0_1 + 100.0*om1_1 + 100.0*omm1_1 - 200.0*om0_1 + 
	100.0*om0_2 + 100.0*om0_0 + 1.5308084989341916e-13*ph1_1 - 
	1.5308084989341916e-13*phm1_1 - 2500.0*ph0_2 + 2500.0*ph0_0)^2 + (-200.0*ps0_1 
	+ 100.0*ps1_1 + 100.0*psm1_1 - 200.0*ps0_1 + 100.0*ps0_2 + 100.0*ps0_0 + 
	0.25*om0_1)^2 + (-25.0*(ps0_2-ps0_0) * (ph1_1-phm1_1) + 25.0*(ps1_1-psm1_1) * 
	(ph0_2-ph0_0) - 200.0*ph0_1 + 100.0*ph1_1 + 100.0*phm1_1 - 200.0*ph0_1 + 
	100.0*ph0_2 + 100.0*ph0_0)^2 + (-200.0*om1_1 + 100.0*om2_1 + 100.0*om0_1 - 
	200.0*om1_1 + 100.0*om1_2 + 100.0*om1_0 + 1.5308084989341916e-13*ph2_1 - 
	1.5308084989341916e-13*ph0_1 - 2500.0*ph1_2 + 2500.0*ph1_0)^2 + (-200.0*ps1_1 + 
	100.0*ps2_1 + 100.0*ps0_1 - 200.0*ps1_1 + 100.0*ps1_2 + 100.0*ps1_0 + 
	0.25*om1_1)^2 + (-25.0*(ps1_2-ps1_0) * (ph2_1-ph0_1) + 25.0*(ps2_1-ps0_1) * 
	(ph1_2-ph1_0) - 200.0*ph1_1 + 100.0*ph2_1 + 100.0*ph0_1 - 200.0*ph1_1 + 
	100.0*ph1_2 + 100.0*ph1_0)^2 + (-200.0*om2_1 + 100.0*om3_1 + 100.0*om1_1 - 
	200.0*om2_1 + 100.0*om2_2 + 100.0*om2_0 + 1.5308084989341916e-13*ph3_1 - 
	1.5308084989341916e-13*ph1_1 - 2500.0*ph2_2 + 2500.0*ph2_0)^2 + (-200.0*ps2_1 + 
	100.0*ps3_1 + 100.0*ps1_1 - 200.0*ps2_1 + 100.0*ps2_2 + 100.0*ps2_0 + 
	0.25*om2_1)^2 + (-25.0*(ps2_2-ps2_0) * (ph3_1-ph1_1) + 25.0*(ps3_1-ps1_1) * 
	(ph2_2-ph2_0) - 200.0*ph2_1 + 100.0*ph3_1 + 100.0*ph1_1 - 200.0*ph2_1 + 
	100.0*ph2_2 + 100.0*ph2_0)^2 + (-200.0*om3_1 + 100.0*om4_1 + 100.0*om2_1 - 
	200.0*om3_1 + 100.0*om3_2 + 100.0*om3_0 + 1.5308084989341916e-13*ph4_1 - 
	1.5308084989341916e-13*ph2_1 - 2500.0*ph3_2 + 2500.0*ph3_0)^2 + (-200.0*ps3_1 + 
	100.0*ps4_1 + 100.0*ps2_1 - 200.0*ps3_1 + 100.0*ps3_2 + 100.0*ps3_0 + 
	0.25*om3_1)^2 + (-25.0*(ps3_2-ps3_0) * (ph4_1-ph2_1) + 25.0*(ps4_1-ps2_1) * 
	(ph3_2-ph3_0) - 200.0*ph3_1 + 100.0*ph4_1 + 100.0*ph2_1 - 200.0*ph3_1 + 
	100.0*ph3_2 + 100.0*ph3_0)^2 + (-200.0*om4_1 + 100.0*om5_1 + 100.0*om3_1 - 
	200.0*om4_1 + 100.0*om4_2 + 100.0*om4_0 + 1.5308084989341916e-13*ph5_1 - 
	1.5308084989341916e-13*ph3_1 - 2500.0*ph4_2 + 2500.0*ph4_0)^2 + (-200.0*ps4_1 + 
	100.0*ps5_1 + 100.0*ps3_1 - 200.0*ps4_1 + 100.0*ps4_2 + 100.0*ps4_0 + 
	0.25*om4_1)^2 + (-25.0*(ps4_2-ps4_0) * (ph5_1-ph3_1) + 25.0*(ps5_1-ps3_1) * 
	(ph4_2-ph4_0) - 200.0*ph4_1 + 100.0*ph5_1 + 100.0*ph3_1 - 200.0*ph4_1 + 
	100.0*ph4_2 + 100.0*ph4_0)^2 + (-200.0*om5_1 + 100.0*om6_1 + 100.0*om4_1 - 
	200.0*om5_1 + 100.0*om5_2 + 100.0*om5_0 + 1.5308084989341916e-13*ph6_1 - 
	1.5308084989341916e-13*ph4_1 - 2500.0*ph5_2 + 2500.0*ph5_0)^2 + (-200.0*ps5_1 + 
	100.0*ps6_1 + 100.0*ps4_1 - 200.0*ps5_1 + 100.0*ps5_2 + 100.0*ps5_0 + 
	0.25*om5_1)^2 + (-25.0*(ps5_2-ps5_0) * (ph6_1-ph4_1) + 25.0*(ps6_1-ps4_1) * 
	(ph5_2-ph5_0) - 200.0*ph5_1 + 100.0*ph6_1 + 100.0*ph4_1 - 200.0*ph5_1 + 
	100.0*ph5_2 + 100.0*ph5_0)^2 + (-200.0*om6_1 + 100.0*om7_1 + 100.0*om5_1 - 
	200.0*om6_1 + 100.0*om6_2 + 100.0*om6_0 + 1.5308084989341916e-13*ph7_1 - 
	1.5308084989341916e-13*ph5_1 - 2500.0*ph6_2 + 2500.0*ph6_0)^2 + (-200.0*ps6_1 + 
	100.0*ps7_1 + 100.0*ps5_1 - 200.0*ps6_1 + 100.0*ps6_2 + 100.0*ps6_0 + 
	0.25*om6_1)^2 + (-25.0*(ps6_2-ps6_0) * (ph7_1-ph5_1) + 25.0*(ps7_1-ps5_1) * 
	(ph6_2-ph6_0) - 200.0*ph6_1 + 100.0*ph7_1 + 100.0*ph5_1 - 200.0*ph6_1 + 
	100.0*ph6_2 + 100.0*ph6_0)^2 + (-200.0*om7_1 + 100.0*om8_1 + 100.0*om6_1 - 
	200.0*om7_1 + 100.0*om7_2 + 100.0*om7_0 + 1.5308084989341916e-13*ph8_1 - 
	1.5308084989341916e-13*ph6_1 - 2500.0*ph7_2 + 2500.0*ph7_0)^2 + (-200.0*ps7_1 + 
	100.0*ps8_1 + 100.0*ps6_1 - 200.0*ps7_1 + 100.0*ps7_2 + 100.0*ps7_0 + 
	0.25*om7_1)^2 + (-25.0*(ps7_2-ps7_0) * (ph8_1-ph6_1) + 25.0*(ps8_1-ps6_1) * 
	(ph7_2-ph7_0) - 200.0*ph7_1 + 100.0*ph8_1 + 100.0*ph6_1 - 200.0*ph7_1 + 
	100.0*ph7_2 + 100.0*ph7_0)^2 + (-200.0*om8_1 + 100.0*om9_1 + 100.0*om7_1 - 
	200.0*om8_1 + 100.0*om8_2 + 100.0*om8_0 + 1.5308084989341916e-13*ph9_1 - 
	1.5308084989341916e-13*ph7_1 - 2500.0*ph8_2 + 2500.0*ph8_0)^2 + (-200.0*ps8_1 + 
	100.0*ps9_1 + 100.0*ps7_1 - 200.0*ps8_1 + 100.0*ps8_2 + 100.0*ps8_0 + 
	0.25*om8_1)^2 + (-25.0*(ps8_2-ps8_0) * (ph9_1-ph7_1) + 25.0*(ps9_1-ps7_1) * 
	(ph8_2-ph8_0) - 200.0*ph8_1 + 100.0*ph9_1 + 100.0*ph7_1 - 200.0*ph8_1 + 
	100.0*ph8_2 + 100.0*ph8_0)^2 + (-200.0*om9_1 + 100.0*om10_1 + 100.0*om8_1 - 
	200.0*om9_1 + 100.0*om9_2 + 100.0*om9_0 + 1.5308084989341916e-13*ph10_1 - 
	1.5308084989341916e-13*ph8_1 - 2500.0*ph9_2 + 2500.0*ph9_0)^2 + (-200.0*ps9_1 + 
	100.0*ps10_1 + 100.0*ps8_1 - 200.0*ps9_1 + 100.0*ps9_2 + 100.0*ps9_0 + 
	0.25*om9_1)^2 + (-25.0*(ps9_2-ps9_0) * (ph10_1-ph8_1) + 25.0*(ps10_1-ps8_1) * 
	(ph9_2-ph9_0) - 200.0*ph9_1 + 100.0*ph10_1 + 100.0*ph8_1 - 200.0*ph9_1 + 
	100.0*ph9_2 + 100.0*ph9_0)^2 + (-200.0*omm9_2 + 100.0*omm8_2 + 100.0*omm10_2 - 
	200.0*omm9_2 + 100.0*omm9_3 + 100.0*omm9_1 + 1.5308084989341916e-13*phm8_2 - 
	1.5308084989341916e-13*phm10_2 - 2500.0*phm9_3 + 2500.0*phm9_1)^2 + 
	(-200.0*psm9_2 + 100.0*psm8_2 + 100.0*psm10_2 - 200.0*psm9_2 + 100.0*psm9_3 + 
	100.0*psm9_1 + 0.25*omm9_2)^2 + (-25.0*(psm9_3-psm9_1) * (phm8_2-phm10_2) + 
	25.0*(psm8_2-psm10_2) * (phm9_3-phm9_1) - 200.0*phm9_2 + 100.0*phm8_2 + 
	100.0*phm10_2 - 200.0*phm9_2 + 100.0*phm9_3 + 100.0*phm9_1)^2 + (-200.0*omm8_2 
	+ 100.0*omm7_2 + 100.0*omm9_2 - 200.0*omm8_2 + 100.0*omm8_3 + 100.0*omm8_1 + 
	1.5308084989341916e-13*phm7_2 - 1.5308084989341916e-13*phm9_2 - 2500.0*phm8_3 + 
	2500.0*phm8_1)^2 + (-200.0*psm8_2 + 100.0*psm7_2 + 100.0*psm9_2 - 200.0*psm8_2 
	+ 100.0*psm8_3 + 100.0*psm8_1 + 0.25*omm8_2)^2 + (-25.0*(psm8_3-psm8_1) * 
	(phm7_2-phm9_2) + 25.0*(psm7_2-psm9_2) * (phm8_3-phm8_1) - 200.0*phm8_2 + 
	100.0*phm7_2 + 100.0*phm9_2 - 200.0*phm8_2 + 100.0*phm8_3 + 100.0*phm8_1)^2 + 
	(-200.0*omm7_2 + 100.0*omm6_2 + 100.0*omm8_2 - 200.0*omm7_2 + 100.0*omm7_3 + 
	100.0*omm7_1 + 1.5308084989341916e-13*phm6_2 - 1.5308084989341916e-13*phm8_2 - 
	2500.0*phm7_3 + 2500.0*phm7_1)^2 + (-200.0*psm7_2 + 100.0*psm6_2 + 100.0*psm8_2 
	- 200.0*psm7_2 + 100.0*psm7_3 + 100.0*psm7_1 + 0.25*omm7_2)^2 + 
	(-25.0*(psm7_3-psm7_1) * (phm6_2-phm8_2) + 25.0*(psm6_2-psm8_2) * 
	(phm7_3-phm7_1) - 200.0*phm7_2 + 100.0*phm6_2 + 100.0*phm8_2 - 200.0*phm7_2 + 
	100.0*phm7_3 + 100.0*phm7_1)^2 + (-200.0*omm6_2 + 100.0*omm5_2 + 100.0*omm7_2 - 
	200.0*omm6_2 + 100.0*omm6_3 + 100.0*omm6_1 + 1.5308084989341916e-13*phm5_2 - 
	1.5308084989341916e-13*phm7_2 - 2500.0*phm6_3 + 2500.0*phm6_1)^2 + 
	(-200.0*psm6_2 + 100.0*psm5_2 + 100.0*psm7_2 - 200.0*psm6_2 + 100.0*psm6_3 + 
	100.0*psm6_1 + 0.25*omm6_2)^2 + (-25.0*(psm6_3-psm6_1) * (phm5_2-phm7_2) + 
	25.0*(psm5_2-psm7_2) * (phm6_3-phm6_1) - 200.0*phm6_2 + 100.0*phm5_2 + 
	100.0*phm7_2 - 200.0*phm6_2 + 100.0*phm6_3 + 100.0*phm6_1)^2 + (-200.0*omm5_2 + 
	100.0*omm4_2 + 100.0*omm6_2 - 200.0*omm5_2 + 100.0*omm5_3 + 100.0*omm5_1 + 
	1.5308084989341916e-13*phm4_2 - 1.5308084989341916e-13*phm6_2 - 2500.0*phm5_3 + 
	2500.0*phm5_1)^2 + (-200.0*psm5_2 + 100.0*psm4_2 + 100.0*psm6_2 - 200.0*psm5_2 
	+ 100.0*psm5_3 + 100.0*psm5_1 + 0.25*omm5_2)^2 + (-25.0*(psm5_3-psm5_1) * 
	(phm4_2-phm6_2) + 25.0*(psm4_2-psm6_2) * (phm5_3-phm5_1) - 200.0*phm5_2 + 
	100.0*phm4_2 + 100.0*phm6_2 - 200.0*phm5_2 + 100.0*phm5_3 + 100.0*phm5_1)^2 + 
	(-200.0*omm4_2 + 100.0*omm3_2 + 100.0*omm5_2 - 200.0*omm4_2 + 100.0*omm4_3 + 
	100.0*omm4_1 + 1.5308084989341916e-13*phm3_2 - 1.5308084989341916e-13*phm5_2 - 
	2500.0*phm4_3 + 2500.0*phm4_1)^2 + (-200.0*psm4_2 + 100.0*psm3_2 + 100.0*psm5_2 
	- 200.0*psm4_2 + 100.0*psm4_3 + 100.0*psm4_1 + 0.25*omm4_2)^2 + 
	(-25.0*(psm4_3-psm4_1) * (phm3_2-phm5_2) + 25.0*(psm3_2-psm5_2) * 
	(phm4_3-phm4_1) - 200.0*phm4_2 + 100.0*phm3_2 + 100.0*phm5_2 - 200.0*phm4_2 + 
	100.0*phm4_3 + 100.0*phm4_1)^2 + (-200.0*omm3_2 + 100.0*omm2_2 + 100.0*omm4_2 - 
	200.0*omm3_2 + 100.0*omm3_3 + 100.0*omm3_1 + 1.5308084989341916e-13*phm2_2 - 
	1.5308084989341916e-13*phm4_2 - 2500.0*phm3_3 + 2500.0*phm3_1)^2 + 
	(-200.0*psm3_2 + 100.0*psm2_2 + 100.0*psm4_2 - 200.0*psm3_2 + 100.0*psm3_3 + 
	100.0*psm3_1 + 0.25*omm3_2)^2 + (-25.0*(psm3_3-psm3_1) * (phm2_2-phm4_2) + 
	25.0*(psm2_2-psm4_2) * (phm3_3-phm3_1) - 200.0*phm3_2 + 100.0*phm2_2 + 
	100.0*phm4_2 - 200.0*phm3_2 + 100.0*phm3_3 + 100.0*phm3_1)^2 + (-200.0*omm2_2 + 
	100.0*omm1_2 + 100.0*omm3_2 - 200.0*omm2_2 + 100.0*omm2_3 + 100.0*omm2_1 + 
	1.5308084989341916e-13*phm1_2 - 1.5308084989341916e-13*phm3_2 - 2500.0*phm2_3 + 
	2500.0*phm2_1)^2 + (-200.0*psm2_2 + 100.0*psm1_2 + 100.0*psm3_2 - 200.0*psm2_2 
	+ 100.0*psm2_3 + 100.0*psm2_1 + 0.25*omm2_2)^2 + (-25.0*(psm2_3-psm2_1) * 
	(phm1_2-phm3_2) + 25.0*(psm1_2-psm3_2) * (phm2_3-phm2_1) - 200.0*phm2_2 + 
	100.0*phm1_2 + 100.0*phm3_2 - 200.0*phm2_2 + 100.0*phm2_3 + 100.0*phm2_1)^2 + 
	(-200.0*omm1_2 + 100.0*om0_2 + 100.0*omm2_2 - 200.0*omm1_2 + 100.0*omm1_3 + 
	100.0*omm1_1 + 1.5308084989341916e-13*ph0_2 - 1.5308084989341916e-13*phm2_2 - 
	2500.0*phm1_3 + 2500.0*phm1_1)^2 + (-200.0*psm1_2 + 100.0*ps0_2 + 100.0*psm2_2 
	- 200.0*psm1_2 + 100.0*psm1_3 + 100.0*psm1_1 + 0.25*omm1_2)^2 + 
	(-25.0*(psm1_3-psm1_1) * (ph0_2-phm2_2) + 25.0*(ps0_2-psm2_2) * (phm1_3-phm1_1) 
	- 200.0*phm1_2 + 100.0*ph0_2 + 100.0*phm2_2 - 200.0*phm1_2 + 100.0*phm1_3 + 
	100.0*phm1_1)^2 + (-200.0*om0_2 + 100.0*om1_2 + 100.0*omm1_2 - 200.0*om0_2 + 
	100.0*om0_3 + 100.0*om0_1 + 1.5308084989341916e-13*ph1_2 - 
	1.5308084989341916e-13*phm1_2 - 2500.0*ph0_3 + 2500.0*ph0_1)^2 + (-200.0*ps0_2 
	+ 100.0*ps1_2 + 100.0*psm1_2 - 200.0*ps0_2 + 100.0*ps0_3 + 100.0*ps0_1 + 
	0.25*om0_2)^2 + (-25.0*(ps0_3-ps0_1) * (ph1_2-phm1_2) + 25.0*(ps1_2-psm1_2) * 
	(ph0_3-ph0_1) - 200.0*ph0_2 + 100.0*ph1_2 + 100.0*phm1_2 - 200.0*ph0_2 + 
	100.0*ph0_3 + 100.0*ph0_1)^2 + (-200.0*om1_2 + 100.0*om2_2 + 100.0*om0_2 - 
	200.0*om1_2 + 100.0*om1_3 + 100.0*om1_1 + 1.5308084989341916e-13*ph2_2 - 
	1.5308084989341916e-13*ph0_2 - 2500.0*ph1_3 + 2500.0*ph1_1)^2 + (-200.0*ps1_2 + 
	100.0*ps2_2 + 100.0*ps0_2 - 200.0*ps1_2 + 100.0*ps1_3 + 100.0*ps1_1 + 
	0.25*om1_2)^2 + (-25.0*(ps1_3-ps1_1) * (ph2_2-ph0_2) + 25.0*(ps2_2-ps0_2) * 
	(ph1_3-ph1_1) - 200.0*ph1_2 + 100.0*ph2_2 + 100.0*ph0_2 - 200.0*ph1_2 + 
	100.0*ph1_3 + 100.0*ph1_1)^2 + (-200.0*om2_2 + 100.0*om3_2 + 100.0*om1_2 - 
	200.0*om2_2 + 100.0*om2_3 + 100.0*om2_1 + 1.5308084989341916e-13*ph3_2 - 
	1.5308084989341916e-13*ph1_2 - 2500.0*ph2_3 + 2500.0*ph2_1)^2 + (-200.0*ps2_2 + 
	100.0*ps3_2 + 100.0*ps1_2 - 200.0*ps2_2 + 100.0*ps2_3 + 100.0*ps2_1 + 
	0.25*om2_2)^2 + (-25.0*(ps2_3-ps2_1) * (ph3_2-ph1_2) + 25.0*(ps3_2-ps1_2) * 
	(ph2_3-ph2_1) - 200.0*ph2_2 + 100.0*ph3_2 + 100.0*ph1_2 - 200.0*ph2_2 + 
	100.0*ph2_3 + 100.0*ph2_1)^2 + (-200.0*om3_2 + 100.0*om4_2 + 100.0*om2_2 - 
	200.0*om3_2 + 100.0*om3_3 + 100.0*om3_1 + 1.5308084989341916e-13*ph4_2 - 
	1.5308084989341916e-13*ph2_2 - 2500.0*ph3_3 + 2500.0*ph3_1)^2 + (-200.0*ps3_2 + 
	100.0*ps4_2 + 100.0*ps2_2 - 200.0*ps3_2 + 100.0*ps3_3 + 100.0*ps3_1 + 
	0.25*om3_2)^2 + (-25.0*(ps3_3-ps3_1) * (ph4_2-ph2_2) + 25.0*(ps4_2-ps2_2) * 
	(ph3_3-ph3_1) - 200.0*ph3_2 + 100.0*ph4_2 + 100.0*ph2_2 - 200.0*ph3_2 + 
	100.0*ph3_3 + 100.0*ph3_1)^2 + (-200.0*om4_2 + 100.0*om5_2 + 100.0*om3_2 - 
	200.0*om4_2 + 100.0*om4_3 + 100.0*om4_1 + 1.5308084989341916e-13*ph5_2 - 
	1.5308084989341916e-13*ph3_2 - 2500.0*ph4_3 + 2500.0*ph4_1)^2 + (-200.0*ps4_2 + 
	100.0*ps5_2 + 100.0*ps3_2 - 200.0*ps4_2 + 100.0*ps4_3 + 100.0*ps4_1 + 
	0.25*om4_2)^2 + (-25.0*(ps4_3-ps4_1) * (ph5_2-ph3_2) + 25.0*(ps5_2-ps3_2) * 
	(ph4_3-ph4_1) - 200.0*ph4_2 + 100.0*ph5_2 + 100.0*ph3_2 - 200.0*ph4_2 + 
	100.0*ph4_3 + 100.0*ph4_1)^2 + (-200.0*om5_2 + 100.0*om6_2 + 100.0*om4_2 - 
	200.0*om5_2 + 100.0*om5_3 + 100.0*om5_1 + 1.5308084989341916e-13*ph6_2 - 
	1.5308084989341916e-13*ph4_2 - 2500.0*ph5_3 + 2500.0*ph5_1)^2 + (-200.0*ps5_2 + 
	100.0*ps6_2 + 100.0*ps4_2 - 200.0*ps5_2 + 100.0*ps5_3 + 100.0*ps5_1 + 
	0.25*om5_2)^2 + (-25.0*(ps5_3-ps5_1) * (ph6_2-ph4_2) + 25.0*(ps6_2-ps4_2) * 
	(ph5_3-ph5_1) - 200.0*ph5_2 + 100.0*ph6_2 + 100.0*ph4_2 - 200.0*ph5_2 + 
	100.0*ph5_3 + 100.0*ph5_1)^2 + (-200.0*om6_2 + 100.0*om7_2 + 100.0*om5_2 - 
	200.0*om6_2 + 100.0*om6_3 + 100.0*om6_1 + 1.5308084989341916e-13*ph7_2 - 
	1.5308084989341916e-13*ph5_2 - 2500.0*ph6_3 + 2500.0*ph6_1)^2 + (-200.0*ps6_2 + 
	100.0*ps7_2 + 100.0*ps5_2 - 200.0*ps6_2 + 100.0*ps6_3 + 100.0*ps6_1 + 
	0.25*om6_2)^2 + (-25.0*(ps6_3-ps6_1) * (ph7_2-ph5_2) + 25.0*(ps7_2-ps5_2) * 
	(ph6_3-ph6_1) - 200.0*ph6_2 + 100.0*ph7_2 + 100.0*ph5_2 - 200.0*ph6_2 + 
	100.0*ph6_3 + 100.0*ph6_1)^2 + (-200.0*om7_2 + 100.0*om8_2 + 100.0*om6_2 - 
	200.0*om7_2 + 100.0*om7_3 + 100.0*om7_1 + 1.5308084989341916e-13*ph8_2 - 
	1.5308084989341916e-13*ph6_2 - 2500.0*ph7_3 + 2500.0*ph7_1)^2 + (-200.0*ps7_2 + 
	100.0*ps8_2 + 100.0*ps6_2 - 200.0*ps7_2 + 100.0*ps7_3 + 100.0*ps7_1 + 
	0.25*om7_2)^2 + (-25.0*(ps7_3-ps7_1) * (ph8_2-ph6_2) + 25.0*(ps8_2-ps6_2) * 
	(ph7_3-ph7_1) - 200.0*ph7_2 + 100.0*ph8_2 + 100.0*ph6_2 - 200.0*ph7_2 + 
	100.0*ph7_3 + 100.0*ph7_1)^2 + (-200.0*om8_2 + 100.0*om9_2 + 100.0*om7_2 - 
	200.0*om8_2 + 100.0*om8_3 + 100.0*om8_1 + 1.5308084989341916e-13*ph9_2 - 
	1.5308084989341916e-13*ph7_2 - 2500.0*ph8_3 + 2500.0*ph8_1)^2 + (-200.0*ps8_2 + 
	100.0*ps9_2 + 100.0*ps7_2 - 200.0*ps8_2 + 100.0*ps8_3 + 100.0*ps8_1 + 
	0.25*om8_2)^2 + (-25.0*(ps8_3-ps8_1) * (ph9_2-ph7_2) + 25.0*(ps9_2-ps7_2) * 
	(ph8_3-ph8_1) - 200.0*ph8_2 + 100.0*ph9_2 + 100.0*ph7_2 - 200.0*ph8_2 + 
	100.0*ph8_3 + 100.0*ph8_1)^2 + (-200.0*om9_2 + 100.0*om10_2 + 100.0*om8_2 - 
	200.0*om9_2 + 100.0*om9_3 + 100.0*om9_1 + 1.5308084989341916e-13*ph10_2 - 
	1.5308084989341916e-13*ph8_2 - 2500.0*ph9_3 + 2500.0*ph9_1)^2 + (-200.0*ps9_2 + 
	100.0*ps10_2 + 100.0*ps8_2 - 200.0*ps9_2 + 100.0*ps9_3 + 100.0*ps9_1 + 
	0.25*om9_2)^2 + (-25.0*(ps9_3-ps9_1) * (ph10_2-ph8_2) + 25.0*(ps10_2-ps8_2) * 
	(ph9_3-ph9_1) - 200.0*ph9_2 + 100.0*ph10_2 + 100.0*ph8_2 - 200.0*ph9_2 + 
	100.0*ph9_3 + 100.0*ph9_1)^2 + (-200.0*omm9_3 + 100.0*omm8_3 + 100.0*omm10_3 - 
	200.0*omm9_3 + 100.0*omm9_4 + 100.0*omm9_2 + 1.5308084989341916e-13*phm8_3 - 
	1.5308084989341916e-13*phm10_3 - 2500.0*phm9_4 + 2500.0*phm9_2)^2 + 
	(-200.0*psm9_3 + 100.0*psm8_3 + 100.0*psm10_3 - 200.0*psm9_3 + 100.0*psm9_4 + 
	100.0*psm9_2 + 0.25*omm9_3)^2 + (-25.0*(psm9_4-psm9_2) * (phm8_3-phm10_3) + 
	25.0*(psm8_3-psm10_3) * (phm9_4-phm9_2) - 200.0*phm9_3 + 100.0*phm8_3 + 
	100.0*phm10_3 - 200.0*phm9_3 + 100.0*phm9_4 + 100.0*phm9_2)^2 + (-200.0*omm8_3 
	+ 100.0*omm7_3 + 100.0*omm9_3 - 200.0*omm8_3 + 100.0*omm8_4 + 100.0*omm8_2 + 
	1.5308084989341916e-13*phm7_3 - 1.5308084989341916e-13*phm9_3 - 2500.0*phm8_4 + 
	2500.0*phm8_2)^2 + (-200.0*psm8_3 + 100.0*psm7_3 + 100.0*psm9_3 - 200.0*psm8_3 
	+ 100.0*psm8_4 + 100.0*psm8_2 + 0.25*omm8_3)^2 + (-25.0*(psm8_4-psm8_2) * 
	(phm7_3-phm9_3) + 25.0*(psm7_3-psm9_3) * (phm8_4-phm8_2) - 200.0*phm8_3 + 
	100.0*phm7_3 + 100.0*phm9_3 - 200.0*phm8_3 + 100.0*phm8_4 + 100.0*phm8_2)^2 + 
	(-200.0*omm7_3 + 100.0*omm6_3 + 100.0*omm8_3 - 200.0*omm7_3 + 100.0*omm7_4 + 
	100.0*omm7_2 + 1.5308084989341916e-13*phm6_3 - 1.5308084989341916e-13*phm8_3 - 
	2500.0*phm7_4 + 2500.0*phm7_2)^2 + (-200.0*psm7_3 + 100.0*psm6_3 + 100.0*psm8_3 
	- 200.0*psm7_3 + 100.0*psm7_4 + 100.0*psm7_2 + 0.25*omm7_3)^2 + 
	(-25.0*(psm7_4-psm7_2) * (phm6_3-phm8_3) + 25.0*(psm6_3-psm8_3) * 
	(phm7_4-phm7_2) - 200.0*phm7_3 + 100.0*phm6_3 + 100.0*phm8_3 - 200.0*phm7_3 + 
	100.0*phm7_4 + 100.0*phm7_2)^2 + (-200.0*omm6_3 + 100.0*omm5_3 + 100.0*omm7_3 - 
	200.0*omm6_3 + 100.0*omm6_4 + 100.0*omm6_2 + 1.5308084989341916e-13*phm5_3 - 
	1.5308084989341916e-13*phm7_3 - 2500.0*phm6_4 + 2500.0*phm6_2)^2 + 
	(-200.0*psm6_3 + 100.0*psm5_3 + 100.0*psm7_3 - 200.0*psm6_3 + 100.0*psm6_4 + 
	100.0*psm6_2 + 0.25*omm6_3)^2 + (-25.0*(psm6_4-psm6_2) * (phm5_3-phm7_3) + 
	25.0*(psm5_3-psm7_3) * (phm6_4-phm6_2) - 200.0*phm6_3 + 100.0*phm5_3 + 
	100.0*phm7_3 - 200.0*phm6_3 + 100.0*phm6_4 + 100.0*phm6_2)^2 + (-200.0*omm5_3 + 
	100.0*omm4_3 + 100.0*omm6_3 - 200.0*omm5_3 + 100.0*omm5_4 + 100.0*omm5_2 + 
	1.5308084989341916e-13*phm4_3 - 1.5308084989341916e-13*phm6_3 - 2500.0*phm5_4 + 
	2500.0*phm5_2)^2 + (-200.0*psm5_3 + 100.0*psm4_3 + 100.0*psm6_3 - 200.0*psm5_3 
	+ 100.0*psm5_4 + 100.0*psm5_2 + 0.25*omm5_3)^2 + (-25.0*(psm5_4-psm5_2) * 
	(phm4_3-phm6_3) + 25.0*(psm4_3-psm6_3) * (phm5_4-phm5_2) - 200.0*phm5_3 + 
	100.0*phm4_3 + 100.0*phm6_3 - 200.0*phm5_3 + 100.0*phm5_4 + 100.0*phm5_2)^2 + 
	(-200.0*omm4_3 + 100.0*omm3_3 + 100.0*omm5_3 - 200.0*omm4_3 + 100.0*omm4_4 + 
	100.0*omm4_2 + 1.5308084989341916e-13*phm3_3 - 1.5308084989341916e-13*phm5_3 - 
	2500.0*phm4_4 + 2500.0*phm4_2)^2 + (-200.0*psm4_3 + 100.0*psm3_3 + 100.0*psm5_3 
	- 200.0*psm4_3 + 100.0*psm4_4 + 100.0*psm4_2 + 0.25*omm4_3)^2 + 
	(-25.0*(psm4_4-psm4_2) * (phm3_3-phm5_3) + 25.0*(psm3_3-psm5_3) * 
	(phm4_4-phm4_2) - 200.0*phm4_3 + 100.0*phm3_3 + 100.0*phm5_3 - 200.0*phm4_3 + 
	100.0*phm4_4 + 100.0*phm4_2)^2 + (-200.0*omm3_3 + 100.0*omm2_3 + 100.0*omm4_3 - 
	200.0*omm3_3 + 100.0*omm3_4 + 100.0*omm3_2 + 1.5308084989341916e-13*phm2_3 - 
	1.5308084989341916e-13*phm4_3 - 2500.0*phm3_4 + 2500.0*phm3_2)^2 + 
	(-200.0*psm3_3 + 100.0*psm2_3 + 100.0*psm4_3 - 200.0*psm3_3 + 100.0*psm3_4 + 
	100.0*psm3_2 + 0.25*omm3_3)^2 + (-25.0*(psm3_4-psm3_2) * (phm2_3-phm4_3) + 
	25.0*(psm2_3-psm4_3) * (phm3_4-phm3_2) - 200.0*phm3_3 + 100.0*phm2_3 + 
	100.0*phm4_3 - 200.0*phm3_3 + 100.0*phm3_4 + 100.0*phm3_2)^2 + (-200.0*omm2_3 + 
	100.0*omm1_3 + 100.0*omm3_3 - 200.0*omm2_3 + 100.0*omm2_4 + 100.0*omm2_2 + 
	1.5308084989341916e-13*phm1_3 - 1.5308084989341916e-13*phm3_3 - 2500.0*phm2_4 + 
	2500.0*phm2_2)^2 + (-200.0*psm2_3 + 100.0*psm1_3 + 100.0*psm3_3 - 200.0*psm2_3 
	+ 100.0*psm2_4 + 100.0*psm2_2 + 0.25*omm2_3)^2 + (-25.0*(psm2_4-psm2_2) * 
	(phm1_3-phm3_3) + 25.0*(psm1_3-psm3_3) * (phm2_4-phm2_2) - 200.0*phm2_3 + 
	100.0*phm1_3 + 100.0*phm3_3 - 200.0*phm2_3 + 100.0*phm2_4 + 100.0*phm2_2)^2 + 
	(-200.0*omm1_3 + 100.0*om0_3 + 100.0*omm2_3 - 200.0*omm1_3 + 100.0*omm1_4 + 
	100.0*omm1_2 + 1.5308084989341916e-13*ph0_3 - 1.5308084989341916e-13*phm2_3 - 
	2500.0*phm1_4 + 2500.0*phm1_2)^2 + (-200.0*psm1_3 + 100.0*ps0_3 + 100.0*psm2_3 
	- 200.0*psm1_3 + 100.0*psm1_4 + 100.0*psm1_2 + 0.25*omm1_3)^2 + 
	(-25.0*(psm1_4-psm1_2) * (ph0_3-phm2_3) + 25.0*(ps0_3-psm2_3) * (phm1_4-phm1_2) 
	- 200.0*phm1_3 + 100.0*ph0_3 + 100.0*phm2_3 - 200.0*phm1_3 + 100.0*phm1_4 + 
	100.0*phm1_2)^2 + (-200.0*om0_3 + 100.0*om1_3 + 100.0*omm1_3 - 200.0*om0_3 + 
	100.0*om0_4 + 100.0*om0_2 + 1.5308084989341916e-13*ph1_3 - 
	1.5308084989341916e-13*phm1_3 - 2500.0*ph0_4 + 2500.0*ph0_2)^2 + (-200.0*ps0_3 
	+ 100.0*ps1_3 + 100.0*psm1_3 - 200.0*ps0_3 + 100.0*ps0_4 + 100.0*ps0_2 + 
	0.25*om0_3)^2 + (-25.0*(ps0_4-ps0_2) * (ph1_3-phm1_3) + 25.0*(ps1_3-psm1_3) * 
	(ph0_4-ph0_2) - 200.0*ph0_3 + 100.0*ph1_3 + 100.0*phm1_3 - 200.0*ph0_3 + 
	100.0*ph0_4 + 100.0*ph0_2)^2 + (-200.0*om1_3 + 100.0*om2_3 + 100.0*om0_3 - 
	200.0*om1_3 + 100.0*om1_4 + 100.0*om1_2 + 1.5308084989341916e-13*ph2_3 - 
	1.5308084989341916e-13*ph0_3 - 2500.0*ph1_4 + 2500.0*ph1_2)^2 + (-200.0*ps1_3 + 
	100.0*ps2_3 + 100.0*ps0_3 - 200.0*ps1_3 + 100.0*ps1_4 + 100.0*ps1_2 + 
	0.25*om1_3)^2 + (-25.0*(ps1_4-ps1_2) * (ph2_3-ph0_3) + 25.0*(ps2_3-ps0_3) * 
	(ph1_4-ph1_2) - 200.0*ph1_3 + 100.0*ph2_3 + 100.0*ph0_3 - 200.0*ph1_3 + 
	100.0*ph1_4 + 100.0*ph1_2)^2 + (-200.0*om2_3 + 100.0*om3_3 + 100.0*om1_3 - 
	200.0*om2_3 + 100.0*om2_4 + 100.0*om2_2 + 1.5308084989341916e-13*ph3_3 - 
	1.5308084989341916e-13*ph1_3 - 2500.0*ph2_4 + 2500.0*ph2_2)^2 + (-200.0*ps2_3 + 
	100.0*ps3_3 + 100.0*ps1_3 - 200.0*ps2_3 + 100.0*ps2_4 + 100.0*ps2_2 + 
	0.25*om2_3)^2 + (-25.0*(ps2_4-ps2_2) * (ph3_3-ph1_3) + 25.0*(ps3_3-ps1_3) * 
	(ph2_4-ph2_2) - 200.0*ph2_3 + 100.0*ph3_3 + 100.0*ph1_3 - 200.0*ph2_3 + 
	100.0*ph2_4 + 100.0*ph2_2)^2 + (-200.0*om3_3 + 100.0*om4_3 + 100.0*om2_3 - 
	200.0*om3_3 + 100.0*om3_4 + 100.0*om3_2 + 1.5308084989341916e-13*ph4_3 - 
	1.5308084989341916e-13*ph2_3 - 2500.0*ph3_4 + 2500.0*ph3_2)^2 + (-200.0*ps3_3 + 
	100.0*ps4_3 + 100.0*ps2_3 - 200.0*ps3_3 + 100.0*ps3_4 + 100.0*ps3_2 + 
	0.25*om3_3)^2 + (-25.0*(ps3_4-ps3_2) * (ph4_3-ph2_3) + 25.0*(ps4_3-ps2_3) * 
	(ph3_4-ph3_2) - 200.0*ph3_3 + 100.0*ph4_3 + 100.0*ph2_3 - 200.0*ph3_3 + 
	100.0*ph3_4 + 100.0*ph3_2)^2 + (-200.0*om4_3 + 100.0*om5_3 + 100.0*om3_3 - 
	200.0*om4_3 + 100.0*om4_4 + 100.0*om4_2 + 1.5308084989341916e-13*ph5_3 - 
	1.5308084989341916e-13*ph3_3 - 2500.0*ph4_4 + 2500.0*ph4_2)^2 + (-200.0*ps4_3 + 
	100.0*ps5_3 + 100.0*ps3_3 - 200.0*ps4_3 + 100.0*ps4_4 + 100.0*ps4_2 + 
	0.25*om4_3)^2 + (-25.0*(ps4_4-ps4_2) * (ph5_3-ph3_3) + 25.0*(ps5_3-ps3_3) * 
	(ph4_4-ph4_2) - 200.0*ph4_3 + 100.0*ph5_3 + 100.0*ph3_3 - 200.0*ph4_3 + 
	100.0*ph4_4 + 100.0*ph4_2)^2 + (-200.0*om5_3 + 100.0*om6_3 + 100.0*om4_3 - 
	200.0*om5_3 + 100.0*om5_4 + 100.0*om5_2 + 1.5308084989341916e-13*ph6_3 - 
	1.5308084989341916e-13*ph4_3 - 2500.0*ph5_4 + 2500.0*ph5_2)^2 + (-200.0*ps5_3 + 
	100.0*ps6_3 + 100.0*ps4_3 - 200.0*ps5_3 + 100.0*ps5_4 + 100.0*ps5_2 + 
	0.25*om5_3)^2 + (-25.0*(ps5_4-ps5_2) * (ph6_3-ph4_3) + 25.0*(ps6_3-ps4_3) * 
	(ph5_4-ph5_2) - 200.0*ph5_3 + 100.0*ph6_3 + 100.0*ph4_3 - 200.0*ph5_3 + 
	100.0*ph5_4 + 100.0*ph5_2)^2 + (-200.0*om6_3 + 100.0*om7_3 + 100.0*om5_3 - 
	200.0*om6_3 + 100.0*om6_4 + 100.0*om6_2 + 1.5308084989341916e-13*ph7_3 - 
	1.5308084989341916e-13*ph5_3 - 2500.0*ph6_4 + 2500.0*ph6_2)^2 + (-200.0*ps6_3 + 
	100.0*ps7_3 + 100.0*ps5_3 - 200.0*ps6_3 + 100.0*ps6_4 + 100.0*ps6_2 + 
	0.25*om6_3)^2 + (-25.0*(ps6_4-ps6_2) * (ph7_3-ph5_3) + 25.0*(ps7_3-ps5_3) * 
	(ph6_4-ph6_2) - 200.0*ph6_3 + 100.0*ph7_3 + 100.0*ph5_3 - 200.0*ph6_3 + 
	100.0*ph6_4 + 100.0*ph6_2)^2 + (-200.0*om7_3 + 100.0*om8_3 + 100.0*om6_3 - 
	200.0*om7_3 + 100.0*om7_4 + 100.0*om7_2 + 1.5308084989341916e-13*ph8_3 - 
	1.5308084989341916e-13*ph6_3 - 2500.0*ph7_4 + 2500.0*ph7_2)^2 + (-200.0*ps7_3 + 
	100.0*ps8_3 + 100.0*ps6_3 - 200.0*ps7_3 + 100.0*ps7_4 + 100.0*ps7_2 + 
	0.25*om7_3)^2 + (-25.0*(ps7_4-ps7_2) * (ph8_3-ph6_3) + 25.0*(ps8_3-ps6_3) * 
	(ph7_4-ph7_2) - 200.0*ph7_3 + 100.0*ph8_3 + 100.0*ph6_3 - 200.0*ph7_3 + 
	100.0*ph7_4 + 100.0*ph7_2)^2 + (-200.0*om8_3 + 100.0*om9_3 + 100.0*om7_3 - 
	200.0*om8_3 + 100.0*om8_4 + 100.0*om8_2 + 1.5308084989341916e-13*ph9_3 - 
	1.5308084989341916e-13*ph7_3 - 2500.0*ph8_4 + 2500.0*ph8_2)^2 + (-200.0*ps8_3 + 
	100.0*ps9_3 + 100.0*ps7_3 - 200.0*ps8_3 + 100.0*ps8_4 + 100.0*ps8_2 + 
	0.25*om8_3)^2 + (-25.0*(ps8_4-ps8_2) * (ph9_3-ph7_3) + 25.0*(ps9_3-ps7_3) * 
	(ph8_4-ph8_2) - 200.0*ph8_3 + 100.0*ph9_3 + 100.0*ph7_3 - 200.0*ph8_3 + 
	100.0*ph8_4 + 100.0*ph8_2)^2 + (-200.0*om9_3 + 100.0*om10_3 + 100.0*om8_3 - 
	200.0*om9_3 + 100.0*om9_4 + 100.0*om9_2 + 1.5308084989341916e-13*ph10_3 - 
	1.5308084989341916e-13*ph8_3 - 2500.0*ph9_4 + 2500.0*ph9_2)^2 + (-200.0*ps9_3 + 
	100.0*ps10_3 + 100.0*ps8_3 - 200.0*ps9_3 + 100.0*ps9_4 + 100.0*ps9_2 + 
	0.25*om9_3)^2 + (-25.0*(ps9_4-ps9_2) * (ph10_3-ph8_3) + 25.0*(ps10_3-ps8_3) * 
	(ph9_4-ph9_2) - 200.0*ph9_3 + 100.0*ph10_3 + 100.0*ph8_3 - 200.0*ph9_3 + 
	100.0*ph9_4 + 100.0*ph9_2)^2 + (-200.0*omm9_4 + 100.0*omm8_4 + 100.0*omm10_4 - 
	200.0*omm9_4 + 100.0*omm9_5 + 100.0*omm9_3 + 1.5308084989341916e-13*phm8_4 - 
	1.5308084989341916e-13*phm10_4 - 2500.0*phm9_5 + 2500.0*phm9_3)^2 + 
	(-200.0*psm9_4 + 100.0*psm8_4 + 100.0*psm10_4 - 200.0*psm9_4 + 100.0*psm9_5 + 
	100.0*psm9_3 + 0.25*omm9_4)^2 + (-25.0*(psm9_5-psm9_3) * (phm8_4-phm10_4) + 
	25.0*(psm8_4-psm10_4) * (phm9_5-phm9_3) - 200.0*phm9_4 + 100.0*phm8_4 + 
	100.0*phm10_4 - 200.0*phm9_4 + 100.0*phm9_5 + 100.0*phm9_3)^2 + (-200.0*omm8_4 
	+ 100.0*omm7_4 + 100.0*omm9_4 - 200.0*omm8_4 + 100.0*omm8_5 + 100.0*omm8_3 + 
	1.5308084989341916e-13*phm7_4 - 1.5308084989341916e-13*phm9_4 - 2500.0*phm8_5 + 
	2500.0*phm8_3)^2 + (-200.0*psm8_4 + 100.0*psm7_4 + 100.0*psm9_4 - 200.0*psm8_4 
	+ 100.0*psm8_5 + 100.0*psm8_3 + 0.25*omm8_4)^2 + (-25.0*(psm8_5-psm8_3) * 
	(phm7_4-phm9_4) + 25.0*(psm7_4-psm9_4) * (phm8_5-phm8_3) - 200.0*phm8_4 + 
	100.0*phm7_4 + 100.0*phm9_4 - 200.0*phm8_4 + 100.0*phm8_5 + 100.0*phm8_3)^2 + 
	(-200.0*omm7_4 + 100.0*omm6_4 + 100.0*omm8_4 - 200.0*omm7_4 + 100.0*omm7_5 + 
	100.0*omm7_3 + 1.5308084989341916e-13*phm6_4 - 1.5308084989341916e-13*phm8_4 - 
	2500.0*phm7_5 + 2500.0*phm7_3)^2 + (-200.0*psm7_4 + 100.0*psm6_4 + 100.0*psm8_4 
	- 200.0*psm7_4 + 100.0*psm7_5 + 100.0*psm7_3 + 0.25*omm7_4)^2 + 
	(-25.0*(psm7_5-psm7_3) * (phm6_4-phm8_4) + 25.0*(psm6_4-psm8_4) * 
	(phm7_5-phm7_3) - 200.0*phm7_4 + 100.0*phm6_4 + 100.0*phm8_4 - 200.0*phm7_4 + 
	100.0*phm7_5 + 100.0*phm7_3)^2 + (-200.0*omm6_4 + 100.0*omm5_4 + 100.0*omm7_4 - 
	200.0*omm6_4 + 100.0*omm6_5 + 100.0*omm6_3 + 1.5308084989341916e-13*phm5_4 - 
	1.5308084989341916e-13*phm7_4 - 2500.0*phm6_5 + 2500.0*phm6_3)^2 + 
	(-200.0*psm6_4 + 100.0*psm5_4 + 100.0*psm7_4 - 200.0*psm6_4 + 100.0*psm6_5 + 
	100.0*psm6_3 + 0.25*omm6_4)^2 + (-25.0*(psm6_5-psm6_3) * (phm5_4-phm7_4) + 
	25.0*(psm5_4-psm7_4) * (phm6_5-phm6_3) - 200.0*phm6_4 + 100.0*phm5_4 + 
	100.0*phm7_4 - 200.0*phm6_4 + 100.0*phm6_5 + 100.0*phm6_3)^2 + (-200.0*omm5_4 + 
	100.0*omm4_4 + 100.0*omm6_4 - 200.0*omm5_4 + 100.0*omm5_5 + 100.0*omm5_3 + 
	1.5308084989341916e-13*phm4_4 - 1.5308084989341916e-13*phm6_4 - 2500.0*phm5_5 + 
	2500.0*phm5_3)^2 + (-200.0*psm5_4 + 100.0*psm4_4 + 100.0*psm6_4 - 200.0*psm5_4 
	+ 100.0*psm5_5 + 100.0*psm5_3 + 0.25*omm5_4)^2 + (-25.0*(psm5_5-psm5_3) * 
	(phm4_4-phm6_4) + 25.0*(psm4_4-psm6_4) * (phm5_5-phm5_3) - 200.0*phm5_4 + 
	100.0*phm4_4 + 100.0*phm6_4 - 200.0*phm5_4 + 100.0*phm5_5 + 100.0*phm5_3)^2 + 
	(-200.0*omm4_4 + 100.0*omm3_4 + 100.0*omm5_4 - 200.0*omm4_4 + 100.0*omm4_5 + 
	100.0*omm4_3 + 1.5308084989341916e-13*phm3_4 - 1.5308084989341916e-13*phm5_4 - 
	2500.0*phm4_5 + 2500.0*phm4_3)^2 + (-200.0*psm4_4 + 100.0*psm3_4 + 100.0*psm5_4 
	- 200.0*psm4_4 + 100.0*psm4_5 + 100.0*psm4_3 + 0.25*omm4_4)^2 + 
	(-25.0*(psm4_5-psm4_3) * (phm3_4-phm5_4) + 25.0*(psm3_4-psm5_4) * 
	(phm4_5-phm4_3) - 200.0*phm4_4 + 100.0*phm3_4 + 100.0*phm5_4 - 200.0*phm4_4 + 
	100.0*phm4_5 + 100.0*phm4_3)^2 + (-200.0*omm3_4 + 100.0*omm2_4 + 100.0*omm4_4 - 
	200.0*omm3_4 + 100.0*omm3_5 + 100.0*omm3_3 + 1.5308084989341916e-13*phm2_4 - 
	1.5308084989341916e-13*phm4_4 - 2500.0*phm3_5 + 2500.0*phm3_3)^2 + 
	(-200.0*psm3_4 + 100.0*psm2_4 + 100.0*psm4_4 - 200.0*psm3_4 + 100.0*psm3_5 + 
	100.0*psm3_3 + 0.25*omm3_4)^2 + (-25.0*(psm3_5-psm3_3) * (phm2_4-phm4_4) + 
	25.0*(psm2_4-psm4_4) * (phm3_5-phm3_3) - 200.0*phm3_4 + 100.0*phm2_4 + 
	100.0*phm4_4 - 200.0*phm3_4 + 100.0*phm3_5 + 100.0*phm3_3)^2 + (-200.0*omm2_4 + 
	100.0*omm1_4 + 100.0*omm3_4 - 200.0*omm2_4 + 100.0*omm2_5 + 100.0*omm2_3 + 
	1.5308084989341916e-13*phm1_4 - 1.5308084989341916e-13*phm3_4 - 2500.0*phm2_5 + 
	2500.0*phm2_3)^2 + (-200.0*psm2_4 + 100.0*psm1_4 + 100.0*psm3_4 - 200.0*psm2_4 
	+ 100.0*psm2_5 + 100.0*psm2_3 + 0.25*omm2_4)^2 + (-25.0*(psm2_5-psm2_3) * 
	(phm1_4-phm3_4) + 25.0*(psm1_4-psm3_4) * (phm2_5-phm2_3) - 200.0*phm2_4 + 
	100.0*phm1_4 + 100.0*phm3_4 - 200.0*phm2_4 + 100.0*phm2_5 + 100.0*phm2_3)^2 + 
	(-200.0*omm1_4 + 100.0*om0_4 + 100.0*omm2_4 - 200.0*omm1_4 + 100.0*omm1_5 + 
	100.0*omm1_3 + 1.5308084989341916e-13*ph0_4 - 1.5308084989341916e-13*phm2_4 - 
	2500.0*phm1_5 + 2500.0*phm1_3)^2 + (-200.0*psm1_4 + 100.0*ps0_4 + 100.0*psm2_4 
	- 200.0*psm1_4 + 100.0*psm1_5 + 100.0*psm1_3 + 0.25*omm1_4)^2 + 
	(-25.0*(psm1_5-psm1_3) * (ph0_4-phm2_4) + 25.0*(ps0_4-psm2_4) * (phm1_5-phm1_3) 
	- 200.0*phm1_4 + 100.0*ph0_4 + 100.0*phm2_4 - 200.0*phm1_4 + 100.0*phm1_5 + 
	100.0*phm1_3)^2 + (-200.0*om0_4 + 100.0*om1_4 + 100.0*omm1_4 - 200.0*om0_4 + 
	100.0*om0_5 + 100.0*om0_3 + 1.5308084989341916e-13*ph1_4 - 
	1.5308084989341916e-13*phm1_4 - 2500.0*ph0_5 + 2500.0*ph0_3)^2 + (-200.0*ps0_4 
	+ 100.0*ps1_4 + 100.0*psm1_4 - 200.0*ps0_4 + 100.0*ps0_5 + 100.0*ps0_3 + 
	0.25*om0_4)^2 + (-25.0*(ps0_5-ps0_3) * (ph1_4-phm1_4) + 25.0*(ps1_4-psm1_4) * 
	(ph0_5-ph0_3) - 200.0*ph0_4 + 100.0*ph1_4 + 100.0*phm1_4 - 200.0*ph0_4 + 
	100.0*ph0_5 + 100.0*ph0_3)^2 + (-200.0*om1_4 + 100.0*om2_4 + 100.0*om0_4 - 
	200.0*om1_4 + 100.0*om1_5 + 100.0*om1_3 + 1.5308084989341916e-13*ph2_4 - 
	1.5308084989341916e-13*ph0_4 - 2500.0*ph1_5 + 2500.0*ph1_3)^2 + (-200.0*ps1_4 + 
	100.0*ps2_4 + 100.0*ps0_4 - 200.0*ps1_4 + 100.0*ps1_5 + 100.0*ps1_3 + 
	0.25*om1_4)^2 + (-25.0*(ps1_5-ps1_3) * (ph2_4-ph0_4) + 25.0*(ps2_4-ps0_4) * 
	(ph1_5-ph1_3) - 200.0*ph1_4 + 100.0*ph2_4 + 100.0*ph0_4 - 200.0*ph1_4 + 
	100.0*ph1_5 + 100.0*ph1_3)^2 + (-200.0*om2_4 + 100.0*om3_4 + 100.0*om1_4 - 
	200.0*om2_4 + 100.0*om2_5 + 100.0*om2_3 + 1.5308084989341916e-13*ph3_4 - 
	1.5308084989341916e-13*ph1_4 - 2500.0*ph2_5 + 2500.0*ph2_3)^2 + (-200.0*ps2_4 + 
	100.0*ps3_4 + 100.0*ps1_4 - 200.0*ps2_4 + 100.0*ps2_5 + 100.0*ps2_3 + 
	0.25*om2_4)^2 + (-25.0*(ps2_5-ps2_3) * (ph3_4-ph1_4) + 25.0*(ps3_4-ps1_4) * 
	(ph2_5-ph2_3) - 200.0*ph2_4 + 100.0*ph3_4 + 100.0*ph1_4 - 200.0*ph2_4 + 
	100.0*ph2_5 + 100.0*ph2_3)^2 + (-200.0*om3_4 + 100.0*om4_4 + 100.0*om2_4 - 
	200.0*om3_4 + 100.0*om3_5 + 100.0*om3_3 + 1.5308084989341916e-13*ph4_4 - 
	1.5308084989341916e-13*ph2_4 - 2500.0*ph3_5 + 2500.0*ph3_3)^2 + (-200.0*ps3_4 + 
	100.0*ps4_4 + 100.0*ps2_4 - 200.0*ps3_4 + 100.0*ps3_5 + 100.0*ps3_3 + 
	0.25*om3_4)^2 + (-25.0*(ps3_5-ps3_3) * (ph4_4-ph2_4) + 25.0*(ps4_4-ps2_4) * 
	(ph3_5-ph3_3) - 200.0*ph3_4 + 100.0*ph4_4 + 100.0*ph2_4 - 200.0*ph3_4 + 
	100.0*ph3_5 + 100.0*ph3_3)^2 + (-200.0*om4_4 + 100.0*om5_4 + 100.0*om3_4 - 
	200.0*om4_4 + 100.0*om4_5 + 100.0*om4_3 + 1.5308084989341916e-13*ph5_4 - 
	1.5308084989341916e-13*ph3_4 - 2500.0*ph4_5 + 2500.0*ph4_3)^2 + (-200.0*ps4_4 + 
	100.0*ps5_4 + 100.0*ps3_4 - 200.0*ps4_4 + 100.0*ps4_5 + 100.0*ps4_3 + 
	0.25*om4_4)^2 + (-25.0*(ps4_5-ps4_3) * (ph5_4-ph3_4) + 25.0*(ps5_4-ps3_4) * 
	(ph4_5-ph4_3) - 200.0*ph4_4 + 100.0*ph5_4 + 100.0*ph3_4 - 200.0*ph4_4 + 
	100.0*ph4_5 + 100.0*ph4_3)^2 + (-200.0*om5_4 + 100.0*om6_4 + 100.0*om4_4 - 
	200.0*om5_4 + 100.0*om5_5 + 100.0*om5_3 + 1.5308084989341916e-13*ph6_4 - 
	1.5308084989341916e-13*ph4_4 - 2500.0*ph5_5 + 2500.0*ph5_3)^2 + (-200.0*ps5_4 + 
	100.0*ps6_4 + 100.0*ps4_4 - 200.0*ps5_4 + 100.0*ps5_5 + 100.0*ps5_3 + 
	0.25*om5_4)^2 + (-25.0*(ps5_5-ps5_3) * (ph6_4-ph4_4) + 25.0*(ps6_4-ps4_4) * 
	(ph5_5-ph5_3) - 200.0*ph5_4 + 100.0*ph6_4 + 100.0*ph4_4 - 200.0*ph5_4 + 
	100.0*ph5_5 + 100.0*ph5_3)^2 + (-200.0*om6_4 + 100.0*om7_4 + 100.0*om5_4 - 
	200.0*om6_4 + 100.0*om6_5 + 100.0*om6_3 + 1.5308084989341916e-13*ph7_4 - 
	1.5308084989341916e-13*ph5_4 - 2500.0*ph6_5 + 2500.0*ph6_3)^2 + (-200.0*ps6_4 + 
	100.0*ps7_4 + 100.0*ps5_4 - 200.0*ps6_4 + 100.0*ps6_5 + 100.0*ps6_3 + 
	0.25*om6_4)^2 + (-25.0*(ps6_5-ps6_3) * (ph7_4-ph5_4) + 25.0*(ps7_4-ps5_4) * 
	(ph6_5-ph6_3) - 200.0*ph6_4 + 100.0*ph7_4 + 100.0*ph5_4 - 200.0*ph6_4 + 
	100.0*ph6_5 + 100.0*ph6_3)^2 + (-200.0*om7_4 + 100.0*om8_4 + 100.0*om6_4 - 
	200.0*om7_4 + 100.0*om7_5 + 100.0*om7_3 + 1.5308084989341916e-13*ph8_4 - 
	1.5308084989341916e-13*ph6_4 - 2500.0*ph7_5 + 2500.0*ph7_3)^2 + (-200.0*ps7_4 + 
	100.0*ps8_4 + 100.0*ps6_4 - 200.0*ps7_4 + 100.0*ps7_5 + 100.0*ps7_3 + 
	0.25*om7_4)^2 + (-25.0*(ps7_5-ps7_3) * (ph8_4-ph6_4) + 25.0*(ps8_4-ps6_4) * 
	(ph7_5-ph7_3) - 200.0*ph7_4 + 100.0*ph8_4 + 100.0*ph6_4 - 200.0*ph7_4 + 
	100.0*ph7_5 + 100.0*ph7_3)^2 + (-200.0*om8_4 + 100.0*om9_4 + 100.0*om7_4 - 
	200.0*om8_4 + 100.0*om8_5 + 100.0*om8_3 + 1.5308084989341916e-13*ph9_4 - 
	1.5308084989341916e-13*ph7_4 - 2500.0*ph8_5 + 2500.0*ph8_3)^2 + (-200.0*ps8_4 + 
	100.0*ps9_4 + 100.0*ps7_4 - 200.0*ps8_4 + 100.0*ps8_5 + 100.0*ps8_3 + 
	0.25*om8_4)^2 + (-25.0*(ps8_5-ps8_3) * (ph9_4-ph7_4) + 25.0*(ps9_4-ps7_4) * 
	(ph8_5-ph8_3) - 200.0*ph8_4 + 100.0*ph9_4 + 100.0*ph7_4 - 200.0*ph8_4 + 
	100.0*ph8_5 + 100.0*ph8_3)^2 + (-200.0*om9_4 + 100.0*om10_4 + 100.0*om8_4 - 
	200.0*om9_4 + 100.0*om9_5 + 100.0*om9_3 + 1.5308084989341916e-13*ph10_4 - 
	1.5308084989341916e-13*ph8_4 - 2500.0*ph9_5 + 2500.0*ph9_3)^2 + (-200.0*ps9_4 + 
	100.0*ps10_4 + 100.0*ps8_4 - 200.0*ps9_4 + 100.0*ps9_5 + 100.0*ps9_3 + 
	0.25*om9_4)^2 + (-25.0*(ps9_5-ps9_3) * (ph10_4-ph8_4) + 25.0*(ps10_4-ps8_4) * 
	(ph9_5-ph9_3) - 200.0*ph9_4 + 100.0*ph10_4 + 100.0*ph8_4 - 200.0*ph9_4 + 
	100.0*ph9_5 + 100.0*ph9_3)^2 + (-200.0*omm9_5 + 100.0*omm8_5 + 100.0*omm10_5 - 
	200.0*omm9_5 + 100.0*omm9_6 + 100.0*omm9_4 + 1.5308084989341916e-13*phm8_5 - 
	1.5308084989341916e-13*phm10_5 - 2500.0*phm9_6 + 2500.0*phm9_4)^2 + 
	(-200.0*psm9_5 + 100.0*psm8_5 + 100.0*psm10_5 - 200.0*psm9_5 + 100.0*psm9_6 + 
	100.0*psm9_4 + 0.25*omm9_5)^2 + (-25.0*(psm9_6-psm9_4) * (phm8_5-phm10_5) + 
	25.0*(psm8_5-psm10_5) * (phm9_6-phm9_4) - 200.0*phm9_5 + 100.0*phm8_5 + 
	100.0*phm10_5 - 200.0*phm9_5 + 100.0*phm9_6 + 100.0*phm9_4)^2 + (-200.0*omm8_5 
	+ 100.0*omm7_5 + 100.0*omm9_5 - 200.0*omm8_5 + 100.0*omm8_6 + 100.0*omm8_4 + 
	1.5308084989341916e-13*phm7_5 - 1.5308084989341916e-13*phm9_5 - 2500.0*phm8_6 + 
	2500.0*phm8_4)^2 + (-200.0*psm8_5 + 100.0*psm7_5 + 100.0*psm9_5 - 200.0*psm8_5 
	+ 100.0*psm8_6 + 100.0*psm8_4 + 0.25*omm8_5)^2 + (-25.0*(psm8_6-psm8_4) * 
	(phm7_5-phm9_5) + 25.0*(psm7_5-psm9_5) * (phm8_6-phm8_4) - 200.0*phm8_5 + 
	100.0*phm7_5 + 100.0*phm9_5 - 200.0*phm8_5 + 100.0*phm8_6 + 100.0*phm8_4)^2 + 
	(-200.0*omm7_5 + 100.0*omm6_5 + 100.0*omm8_5 - 200.0*omm7_5 + 100.0*omm7_6 + 
	100.0*omm7_4 + 1.5308084989341916e-13*phm6_5 - 1.5308084989341916e-13*phm8_5 - 
	2500.0*phm7_6 + 2500.0*phm7_4)^2 + (-200.0*psm7_5 + 100.0*psm6_5 + 100.0*psm8_5 
	- 200.0*psm7_5 + 100.0*psm7_6 + 100.0*psm7_4 + 0.25*omm7_5)^2 + 
	(-25.0*(psm7_6-psm7_4) * (phm6_5-phm8_5) + 25.0*(psm6_5-psm8_5) * 
	(phm7_6-phm7_4) - 200.0*phm7_5 + 100.0*phm6_5 + 100.0*phm8_5 - 200.0*phm7_5 + 
	100.0*phm7_6 + 100.0*phm7_4)^2 + (-200.0*omm6_5 + 100.0*omm5_5 + 100.0*omm7_5 - 
	200.0*omm6_5 + 100.0*omm6_6 + 100.0*omm6_4 + 1.5308084989341916e-13*phm5_5 - 
	1.5308084989341916e-13*phm7_5 - 2500.0*phm6_6 + 2500.0*phm6_4)^2 + 
	(-200.0*psm6_5 + 100.0*psm5_5 + 100.0*psm7_5 - 200.0*psm6_5 + 100.0*psm6_6 + 
	100.0*psm6_4 + 0.25*omm6_5)^2 + (-25.0*(psm6_6-psm6_4) * (phm5_5-phm7_5) + 
	25.0*(psm5_5-psm7_5) * (phm6_6-phm6_4) - 200.0*phm6_5 + 100.0*phm5_5 + 
	100.0*phm7_5 - 200.0*phm6_5 + 100.0*phm6_6 + 100.0*phm6_4)^2 + (-200.0*omm5_5 + 
	100.0*omm4_5 + 100.0*omm6_5 - 200.0*omm5_5 + 100.0*omm5_6 + 100.0*omm5_4 + 
	1.5308084989341916e-13*phm4_5 - 1.5308084989341916e-13*phm6_5 - 2500.0*phm5_6 + 
	2500.0*phm5_4)^2 + (-200.0*psm5_5 + 100.0*psm4_5 + 100.0*psm6_5 - 200.0*psm5_5 
	+ 100.0*psm5_6 + 100.0*psm5_4 + 0.25*omm5_5)^2 + (-25.0*(psm5_6-psm5_4) * 
	(phm4_5-phm6_5) + 25.0*(psm4_5-psm6_5) * (phm5_6-phm5_4) - 200.0*phm5_5 + 
	100.0*phm4_5 + 100.0*phm6_5 - 200.0*phm5_5 + 100.0*phm5_6 + 100.0*phm5_4)^2 + 
	(-200.0*omm4_5 + 100.0*omm3_5 + 100.0*omm5_5 - 200.0*omm4_5 + 100.0*omm4_6 + 
	100.0*omm4_4 + 1.5308084989341916e-13*phm3_5 - 1.5308084989341916e-13*phm5_5 - 
	2500.0*phm4_6 + 2500.0*phm4_4)^2 + (-200.0*psm4_5 + 100.0*psm3_5 + 100.0*psm5_5 
	- 200.0*psm4_5 + 100.0*psm4_6 + 100.0*psm4_4 + 0.25*omm4_5)^2 + 
	(-25.0*(psm4_6-psm4_4) * (phm3_5-phm5_5) + 25.0*(psm3_5-psm5_5) * 
	(phm4_6-phm4_4) - 200.0*phm4_5 + 100.0*phm3_5 + 100.0*phm5_5 - 200.0*phm4_5 + 
	100.0*phm4_6 + 100.0*phm4_4)^2 + (-200.0*omm3_5 + 100.0*omm2_5 + 100.0*omm4_5 - 
	200.0*omm3_5 + 100.0*omm3_6 + 100.0*omm3_4 + 1.5308084989341916e-13*phm2_5 - 
	1.5308084989341916e-13*phm4_5 - 2500.0*phm3_6 + 2500.0*phm3_4)^2 + 
	(-200.0*psm3_5 + 100.0*psm2_5 + 100.0*psm4_5 - 200.0*psm3_5 + 100.0*psm3_6 + 
	100.0*psm3_4 + 0.25*omm3_5)^2 + (-25.0*(psm3_6-psm3_4) * (phm2_5-phm4_5) + 
	25.0*(psm2_5-psm4_5) * (phm3_6-phm3_4) - 200.0*phm3_5 + 100.0*phm2_5 + 
	100.0*phm4_5 - 200.0*phm3_5 + 100.0*phm3_6 + 100.0*phm3_4)^2 + (-200.0*omm2_5 + 
	100.0*omm1_5 + 100.0*omm3_5 - 200.0*omm2_5 + 100.0*omm2_6 + 100.0*omm2_4 + 
	1.5308084989341916e-13*phm1_5 - 1.5308084989341916e-13*phm3_5 - 2500.0*phm2_6 + 
	2500.0*phm2_4)^2 + (-200.0*psm2_5 + 100.0*psm1_5 + 100.0*psm3_5 - 200.0*psm2_5 
	+ 100.0*psm2_6 + 100.0*psm2_4 + 0.25*omm2_5)^2 + (-25.0*(psm2_6-psm2_4) * 
	(phm1_5-phm3_5) + 25.0*(psm1_5-psm3_5) * (phm2_6-phm2_4) - 200.0*phm2_5 + 
	100.0*phm1_5 + 100.0*phm3_5 - 200.0*phm2_5 + 100.0*phm2_6 + 100.0*phm2_4)^2 + 
	(-200.0*omm1_5 + 100.0*om0_5 + 100.0*omm2_5 - 200.0*omm1_5 + 100.0*omm1_6 + 
	100.0*omm1_4 + 1.5308084989341916e-13*ph0_5 - 1.5308084989341916e-13*phm2_5 - 
	2500.0*phm1_6 + 2500.0*phm1_4)^2 + (-200.0*psm1_5 + 100.0*ps0_5 + 100.0*psm2_5 
	- 200.0*psm1_5 + 100.0*psm1_6 + 100.0*psm1_4 + 0.25*omm1_5)^2 + 
	(-25.0*(psm1_6-psm1_4) * (ph0_5-phm2_5) + 25.0*(ps0_5-psm2_5) * (phm1_6-phm1_4) 
	- 200.0*phm1_5 + 100.0*ph0_5 + 100.0*phm2_5 - 200.0*phm1_5 + 100.0*phm1_6 + 
	100.0*phm1_4)^2 + (-200.0*om0_5 + 100.0*om1_5 + 100.0*omm1_5 - 200.0*om0_5 + 
	100.0*om0_6 + 100.0*om0_4 + 1.5308084989341916e-13*ph1_5 - 
	1.5308084989341916e-13*phm1_5 - 2500.0*ph0_6 + 2500.0*ph0_4)^2 + (-200.0*ps0_5 
	+ 100.0*ps1_5 + 100.0*psm1_5 - 200.0*ps0_5 + 100.0*ps0_6 + 100.0*ps0_4 + 
	0.25*om0_5)^2 + (-25.0*(ps0_6-ps0_4) * (ph1_5-phm1_5) + 25.0*(ps1_5-psm1_5) * 
	(ph0_6-ph0_4) - 200.0*ph0_5 + 100.0*ph1_5 + 100.0*phm1_5 - 200.0*ph0_5 + 
	100.0*ph0_6 + 100.0*ph0_4)^2 + (-200.0*om1_5 + 100.0*om2_5 + 100.0*om0_5 - 
	200.0*om1_5 + 100.0*om1_6 + 100.0*om1_4 + 1.5308084989341916e-13*ph2_5 - 
	1.5308084989341916e-13*ph0_5 - 2500.0*ph1_6 + 2500.0*ph1_4)^2 + (-200.0*ps1_5 + 
	100.0*ps2_5 + 100.0*ps0_5 - 200.0*ps1_5 + 100.0*ps1_6 + 100.0*ps1_4 + 
	0.25*om1_5)^2 + (-25.0*(ps1_6-ps1_4) * (ph2_5-ph0_5) + 25.0*(ps2_5-ps0_5) * 
	(ph1_6-ph1_4) - 200.0*ph1_5 + 100.0*ph2_5 + 100.0*ph0_5 - 200.0*ph1_5 + 
	100.0*ph1_6 + 100.0*ph1_4)^2 + (-200.0*om2_5 + 100.0*om3_5 + 100.0*om1_5 - 
	200.0*om2_5 + 100.0*om2_6 + 100.0*om2_4 + 1.5308084989341916e-13*ph3_5 - 
	1.5308084989341916e-13*ph1_5 - 2500.0*ph2_6 + 2500.0*ph2_4)^2 + (-200.0*ps2_5 + 
	100.0*ps3_5 + 100.0*ps1_5 - 200.0*ps2_5 + 100.0*ps2_6 + 100.0*ps2_4 + 
	0.25*om2_5)^2 + (-25.0*(ps2_6-ps2_4) * (ph3_5-ph1_5) + 25.0*(ps3_5-ps1_5) * 
	(ph2_6-ph2_4) - 200.0*ph2_5 + 100.0*ph3_5 + 100.0*ph1_5 - 200.0*ph2_5 + 
	100.0*ph2_6 + 100.0*ph2_4)^2 + (-200.0*om3_5 + 100.0*om4_5 + 100.0*om2_5 - 
	200.0*om3_5 + 100.0*om3_6 + 100.0*om3_4 + 1.5308084989341916e-13*ph4_5 - 
	1.5308084989341916e-13*ph2_5 - 2500.0*ph3_6 + 2500.0*ph3_4)^2 + (-200.0*ps3_5 + 
	100.0*ps4_5 + 100.0*ps2_5 - 200.0*ps3_5 + 100.0*ps3_6 + 100.0*ps3_4 + 
	0.25*om3_5)^2 + (-25.0*(ps3_6-ps3_4) * (ph4_5-ph2_5) + 25.0*(ps4_5-ps2_5) * 
	(ph3_6-ph3_4) - 200.0*ph3_5 + 100.0*ph4_5 + 100.0*ph2_5 - 200.0*ph3_5 + 
	100.0*ph3_6 + 100.0*ph3_4)^2 + (-200.0*om4_5 + 100.0*om5_5 + 100.0*om3_5 - 
	200.0*om4_5 + 100.0*om4_6 + 100.0*om4_4 + 1.5308084989341916e-13*ph5_5 - 
	1.5308084989341916e-13*ph3_5 - 2500.0*ph4_6 + 2500.0*ph4_4)^2 + (-200.0*ps4_5 + 
	100.0*ps5_5 + 100.0*ps3_5 - 200.0*ps4_5 + 100.0*ps4_6 + 100.0*ps4_4 + 
	0.25*om4_5)^2 + (-25.0*(ps4_6-ps4_4) * (ph5_5-ph3_5) + 25.0*(ps5_5-ps3_5) * 
	(ph4_6-ph4_4) - 200.0*ph4_5 + 100.0*ph5_5 + 100.0*ph3_5 - 200.0*ph4_5 + 
	100.0*ph4_6 + 100.0*ph4_4)^2 + (-200.0*om5_5 + 100.0*om6_5 + 100.0*om4_5 - 
	200.0*om5_5 + 100.0*om5_6 + 100.0*om5_4 + 1.5308084989341916e-13*ph6_5 - 
	1.5308084989341916e-13*ph4_5 - 2500.0*ph5_6 + 2500.0*ph5_4)^2 + (-200.0*ps5_5 + 
	100.0*ps6_5 + 100.0*ps4_5 - 200.0*ps5_5 + 100.0*ps5_6 + 100.0*ps5_4 + 
	0.25*om5_5)^2 + (-25.0*(ps5_6-ps5_4) * (ph6_5-ph4_5) + 25.0*(ps6_5-ps4_5) * 
	(ph5_6-ph5_4) - 200.0*ph5_5 + 100.0*ph6_5 + 100.0*ph4_5 - 200.0*ph5_5 + 
	100.0*ph5_6 + 100.0*ph5_4)^2 + (-200.0*om6_5 + 100.0*om7_5 + 100.0*om5_5 - 
	200.0*om6_5 + 100.0*om6_6 + 100.0*om6_4 + 1.5308084989341916e-13*ph7_5 - 
	1.5308084989341916e-13*ph5_5 - 2500.0*ph6_6 + 2500.0*ph6_4)^2 + (-200.0*ps6_5 + 
	100.0*ps7_5 + 100.0*ps5_5 - 200.0*ps6_5 + 100.0*ps6_6 + 100.0*ps6_4 + 
	0.25*om6_5)^2 + (-25.0*(ps6_6-ps6_4) * (ph7_5-ph5_5) + 25.0*(ps7_5-ps5_5) * 
	(ph6_6-ph6_4) - 200.0*ph6_5 + 100.0*ph7_5 + 100.0*ph5_5 - 200.0*ph6_5 + 
	100.0*ph6_6 + 100.0*ph6_4)^2 + (-200.0*om7_5 + 100.0*om8_5 + 100.0*om6_5 - 
	200.0*om7_5 + 100.0*om7_6 + 100.0*om7_4 + 1.5308084989341916e-13*ph8_5 - 
	1.5308084989341916e-13*ph6_5 - 2500.0*ph7_6 + 2500.0*ph7_4)^2 + (-200.0*ps7_5 + 
	100.0*ps8_5 + 100.0*ps6_5 - 200.0*ps7_5 + 100.0*ps7_6 + 100.0*ps7_4 + 
	0.25*om7_5)^2 + (-25.0*(ps7_6-ps7_4) * (ph8_5-ph6_5) + 25.0*(ps8_5-ps6_5) * 
	(ph7_6-ph7_4) - 200.0*ph7_5 + 100.0*ph8_5 + 100.0*ph6_5 - 200.0*ph7_5 + 
	100.0*ph7_6 + 100.0*ph7_4)^2 + (-200.0*om8_5 + 100.0*om9_5 + 100.0*om7_5 - 
	200.0*om8_5 + 100.0*om8_6 + 100.0*om8_4 + 1.5308084989341916e-13*ph9_5 - 
	1.5308084989341916e-13*ph7_5 - 2500.0*ph8_6 + 2500.0*ph8_4)^2 + (-200.0*ps8_5 + 
	100.0*ps9_5 + 100.0*ps7_5 - 200.0*ps8_5 + 100.0*ps8_6 + 100.0*ps8_4 + 
	0.25*om8_5)^2 + (-25.0*(ps8_6-ps8_4) * (ph9_5-ph7_5) + 25.0*(ps9_5-ps7_5) * 
	(ph8_6-ph8_4) - 200.0*ph8_5 + 100.0*ph9_5 + 100.0*ph7_5 - 200.0*ph8_5 + 
	100.0*ph8_6 + 100.0*ph8_4)^2 + (-200.0*om9_5 + 100.0*om10_5 + 100.0*om8_5 - 
	200.0*om9_5 + 100.0*om9_6 + 100.0*om9_4 + 1.5308084989341916e-13*ph10_5 - 
	1.5308084989341916e-13*ph8_5 - 2500.0*ph9_6 + 2500.0*ph9_4)^2 + (-200.0*ps9_5 + 
	100.0*ps10_5 + 100.0*ps8_5 - 200.0*ps9_5 + 100.0*ps9_6 + 100.0*ps9_4 + 
	0.25*om9_5)^2 + (-25.0*(ps9_6-ps9_4) * (ph10_5-ph8_5) + 25.0*(ps10_5-ps8_5) * 
	(ph9_6-ph9_4) - 200.0*ph9_5 + 100.0*ph10_5 + 100.0*ph8_5 - 200.0*ph9_5 + 
	100.0*ph9_6 + 100.0*ph9_4)^2 + (-200.0*omm9_6 + 100.0*omm8_6 + 100.0*omm10_6 - 
	200.0*omm9_6 + 100.0*omm9_7 + 100.0*omm9_5 + 1.5308084989341916e-13*phm8_6 - 
	1.5308084989341916e-13*phm10_6 - 2500.0*phm9_7 + 2500.0*phm9_5)^2 + 
	(-200.0*psm9_6 + 100.0*psm8_6 + 100.0*psm10_6 - 200.0*psm9_6 + 100.0*psm9_7 + 
	100.0*psm9_5 + 0.25*omm9_6)^2 + (-25.0*(psm9_7-psm9_5) * (phm8_6-phm10_6) + 
	25.0*(psm8_6-psm10_6) * (phm9_7-phm9_5) - 200.0*phm9_6 + 100.0*phm8_6 + 
	100.0*phm10_6 - 200.0*phm9_6 + 100.0*phm9_7 + 100.0*phm9_5)^2 + (-200.0*omm8_6 
	+ 100.0*omm7_6 + 100.0*omm9_6 - 200.0*omm8_6 + 100.0*omm8_7 + 100.0*omm8_5 + 
	1.5308084989341916e-13*phm7_6 - 1.5308084989341916e-13*phm9_6 - 2500.0*phm8_7 + 
	2500.0*phm8_5)^2 + (-200.0*psm8_6 + 100.0*psm7_6 + 100.0*psm9_6 - 200.0*psm8_6 
	+ 100.0*psm8_7 + 100.0*psm8_5 + 0.25*omm8_6)^2 + (-25.0*(psm8_7-psm8_5) * 
	(phm7_6-phm9_6) + 25.0*(psm7_6-psm9_6) * (phm8_7-phm8_5) - 200.0*phm8_6 + 
	100.0*phm7_6 + 100.0*phm9_6 - 200.0*phm8_6 + 100.0*phm8_7 + 100.0*phm8_5)^2 + 
	(-200.0*omm7_6 + 100.0*omm6_6 + 100.0*omm8_6 - 200.0*omm7_6 + 100.0*omm7_7 + 
	100.0*omm7_5 + 1.5308084989341916e-13*phm6_6 - 1.5308084989341916e-13*phm8_6 - 
	2500.0*phm7_7 + 2500.0*phm7_5)^2 + (-200.0*psm7_6 + 100.0*psm6_6 + 100.0*psm8_6 
	- 200.0*psm7_6 + 100.0*psm7_7 + 100.0*psm7_5 + 0.25*omm7_6)^2 + 
	(-25.0*(psm7_7-psm7_5) * (phm6_6-phm8_6) + 25.0*(psm6_6-psm8_6) * 
	(phm7_7-phm7_5) - 200.0*phm7_6 + 100.0*phm6_6 + 100.0*phm8_6 - 200.0*phm7_6 + 
	100.0*phm7_7 + 100.0*phm7_5)^2 + (-200.0*omm6_6 + 100.0*omm5_6 + 100.0*omm7_6 - 
	200.0*omm6_6 + 100.0*omm6_7 + 100.0*omm6_5 + 1.5308084989341916e-13*phm5_6 - 
	1.5308084989341916e-13*phm7_6 - 2500.0*phm6_7 + 2500.0*phm6_5)^2 + 
	(-200.0*psm6_6 + 100.0*psm5_6 + 100.0*psm7_6 - 200.0*psm6_6 + 100.0*psm6_7 + 
	100.0*psm6_5 + 0.25*omm6_6)^2 + (-25.0*(psm6_7-psm6_5) * (phm5_6-phm7_6) + 
	25.0*(psm5_6-psm7_6) * (phm6_7-phm6_5) - 200.0*phm6_6 + 100.0*phm5_6 + 
	100.0*phm7_6 - 200.0*phm6_6 + 100.0*phm6_7 + 100.0*phm6_5)^2 + (-200.0*omm5_6 + 
	100.0*omm4_6 + 100.0*omm6_6 - 200.0*omm5_6 + 100.0*omm5_7 + 100.0*omm5_5 + 
	1.5308084989341916e-13*phm4_6 - 1.5308084989341916e-13*phm6_6 - 2500.0*phm5_7 + 
	2500.0*phm5_5)^2 + (-200.0*psm5_6 + 100.0*psm4_6 + 100.0*psm6_6 - 200.0*psm5_6 
	+ 100.0*psm5_7 + 100.0*psm5_5 + 0.25*omm5_6)^2 + (-25.0*(psm5_7-psm5_5) * 
	(phm4_6-phm6_6) + 25.0*(psm4_6-psm6_6) * (phm5_7-phm5_5) - 200.0*phm5_6 + 
	100.0*phm4_6 + 100.0*phm6_6 - 200.0*phm5_6 + 100.0*phm5_7 + 100.0*phm5_5)^2 + 
	(-200.0*omm4_6 + 100.0*omm3_6 + 100.0*omm5_6 - 200.0*omm4_6 + 100.0*omm4_7 + 
	100.0*omm4_5 + 1.5308084989341916e-13*phm3_6 - 1.5308084989341916e-13*phm5_6 - 
	2500.0*phm4_7 + 2500.0*phm4_5)^2 + (-200.0*psm4_6 + 100.0*psm3_6 + 100.0*psm5_6 
	- 200.0*psm4_6 + 100.0*psm4_7 + 100.0*psm4_5 + 0.25*omm4_6)^2 + 
	(-25.0*(psm4_7-psm4_5) * (phm3_6-phm5_6) + 25.0*(psm3_6-psm5_6) * 
	(phm4_7-phm4_5) - 200.0*phm4_6 + 100.0*phm3_6 + 100.0*phm5_6 - 200.0*phm4_6 + 
	100.0*phm4_7 + 100.0*phm4_5)^2 + (-200.0*omm3_6 + 100.0*omm2_6 + 100.0*omm4_6 - 
	200.0*omm3_6 + 100.0*omm3_7 + 100.0*omm3_5 + 1.5308084989341916e-13*phm2_6 - 
	1.5308084989341916e-13*phm4_6 - 2500.0*phm3_7 + 2500.0*phm3_5)^2 + 
	(-200.0*psm3_6 + 100.0*psm2_6 + 100.0*psm4_6 - 200.0*psm3_6 + 100.0*psm3_7 + 
	100.0*psm3_5 + 0.25*omm3_6)^2 + (-25.0*(psm3_7-psm3_5) * (phm2_6-phm4_6) + 
	25.0*(psm2_6-psm4_6) * (phm3_7-phm3_5) - 200.0*phm3_6 + 100.0*phm2_6 + 
	100.0*phm4_6 - 200.0*phm3_6 + 100.0*phm3_7 + 100.0*phm3_5)^2 + (-200.0*omm2_6 + 
	100.0*omm1_6 + 100.0*omm3_6 - 200.0*omm2_6 + 100.0*omm2_7 + 100.0*omm2_5 + 
	1.5308084989341916e-13*phm1_6 - 1.5308084989341916e-13*phm3_6 - 2500.0*phm2_7 + 
	2500.0*phm2_5)^2 + (-200.0*psm2_6 + 100.0*psm1_6 + 100.0*psm3_6 - 200.0*psm2_6 
	+ 100.0*psm2_7 + 100.0*psm2_5 + 0.25*omm2_6)^2 + (-25.0*(psm2_7-psm2_5) * 
	(phm1_6-phm3_6) + 25.0*(psm1_6-psm3_6) * (phm2_7-phm2_5) - 200.0*phm2_6 + 
	100.0*phm1_6 + 100.0*phm3_6 - 200.0*phm2_6 + 100.0*phm2_7 + 100.0*phm2_5)^2 + 
	(-200.0*omm1_6 + 100.0*om0_6 + 100.0*omm2_6 - 200.0*omm1_6 + 100.0*omm1_7 + 
	100.0*omm1_5 + 1.5308084989341916e-13*ph0_6 - 1.5308084989341916e-13*phm2_6 - 
	2500.0*phm1_7 + 2500.0*phm1_5)^2 + (-200.0*psm1_6 + 100.0*ps0_6 + 100.0*psm2_6 
	- 200.0*psm1_6 + 100.0*psm1_7 + 100.0*psm1_5 + 0.25*omm1_6)^2 + 
	(-25.0*(psm1_7-psm1_5) * (ph0_6-phm2_6) + 25.0*(ps0_6-psm2_6) * (phm1_7-phm1_5) 
	- 200.0*phm1_6 + 100.0*ph0_6 + 100.0*phm2_6 - 200.0*phm1_6 + 100.0*phm1_7 + 
	100.0*phm1_5)^2 + (-200.0*om0_6 + 100.0*om1_6 + 100.0*omm1_6 - 200.0*om0_6 + 
	100.0*om0_7 + 100.0*om0_5 + 1.5308084989341916e-13*ph1_6 - 
	1.5308084989341916e-13*phm1_6 - 2500.0*ph0_7 + 2500.0*ph0_5)^2 + (-200.0*ps0_6 
	+ 100.0*ps1_6 + 100.0*psm1_6 - 200.0*ps0_6 + 100.0*ps0_7 + 100.0*ps0_5 + 
	0.25*om0_6)^2 + (-25.0*(ps0_7-ps0_5) * (ph1_6-phm1_6) + 25.0*(ps1_6-psm1_6) * 
	(ph0_7-ph0_5) - 200.0*ph0_6 + 100.0*ph1_6 + 100.0*phm1_6 - 200.0*ph0_6 + 
	100.0*ph0_7 + 100.0*ph0_5)^2 + (-200.0*om1_6 + 100.0*om2_6 + 100.0*om0_6 - 
	200.0*om1_6 + 100.0*om1_7 + 100.0*om1_5 + 1.5308084989341916e-13*ph2_6 - 
	1.5308084989341916e-13*ph0_6 - 2500.0*ph1_7 + 2500.0*ph1_5)^2 + (-200.0*ps1_6 + 
	100.0*ps2_6 + 100.0*ps0_6 - 200.0*ps1_6 + 100.0*ps1_7 + 100.0*ps1_5 + 
	0.25*om1_6)^2 + (-25.0*(ps1_7-ps1_5) * (ph2_6-ph0_6) + 25.0*(ps2_6-ps0_6) * 
	(ph1_7-ph1_5) - 200.0*ph1_6 + 100.0*ph2_6 + 100.0*ph0_6 - 200.0*ph1_6 + 
	100.0*ph1_7 + 100.0*ph1_5)^2 + (-200.0*om2_6 + 100.0*om3_6 + 100.0*om1_6 - 
	200.0*om2_6 + 100.0*om2_7 + 100.0*om2_5 + 1.5308084989341916e-13*ph3_6 - 
	1.5308084989341916e-13*ph1_6 - 2500.0*ph2_7 + 2500.0*ph2_5)^2 + (-200.0*ps2_6 + 
	100.0*ps3_6 + 100.0*ps1_6 - 200.0*ps2_6 + 100.0*ps2_7 + 100.0*ps2_5 + 
	0.25*om2_6)^2 + (-25.0*(ps2_7-ps2_5) * (ph3_6-ph1_6) + 25.0*(ps3_6-ps1_6) * 
	(ph2_7-ph2_5) - 200.0*ph2_6 + 100.0*ph3_6 + 100.0*ph1_6 - 200.0*ph2_6 + 
	100.0*ph2_7 + 100.0*ph2_5)^2 + (-200.0*om3_6 + 100.0*om4_6 + 100.0*om2_6 - 
	200.0*om3_6 + 100.0*om3_7 + 100.0*om3_5 + 1.5308084989341916e-13*ph4_6 - 
	1.5308084989341916e-13*ph2_6 - 2500.0*ph3_7 + 2500.0*ph3_5)^2 + (-200.0*ps3_6 + 
	100.0*ps4_6 + 100.0*ps2_6 - 200.0*ps3_6 + 100.0*ps3_7 + 100.0*ps3_5 + 
	0.25*om3_6)^2 + (-25.0*(ps3_7-ps3_5) * (ph4_6-ph2_6) + 25.0*(ps4_6-ps2_6) * 
	(ph3_7-ph3_5) - 200.0*ph3_6 + 100.0*ph4_6 + 100.0*ph2_6 - 200.0*ph3_6 + 
	100.0*ph3_7 + 100.0*ph3_5)^2 + (-200.0*om4_6 + 100.0*om5_6 + 100.0*om3_6 - 
	200.0*om4_6 + 100.0*om4_7 + 100.0*om4_5 + 1.5308084989341916e-13*ph5_6 - 
	1.5308084989341916e-13*ph3_6 - 2500.0*ph4_7 + 2500.0*ph4_5)^2 + (-200.0*ps4_6 + 
	100.0*ps5_6 + 100.0*ps3_6 - 200.0*ps4_6 + 100.0*ps4_7 + 100.0*ps4_5 + 
	0.25*om4_6)^2 + (-25.0*(ps4_7-ps4_5) * (ph5_6-ph3_6) + 25.0*(ps5_6-ps3_6) * 
	(ph4_7-ph4_5) - 200.0*ph4_6 + 100.0*ph5_6 + 100.0*ph3_6 - 200.0*ph4_6 + 
	100.0*ph4_7 + 100.0*ph4_5)^2 + (-200.0*om5_6 + 100.0*om6_6 + 100.0*om4_6 - 
	200.0*om5_6 + 100.0*om5_7 + 100.0*om5_5 + 1.5308084989341916e-13*ph6_6 - 
	1.5308084989341916e-13*ph4_6 - 2500.0*ph5_7 + 2500.0*ph5_5)^2 + (-200.0*ps5_6 + 
	100.0*ps6_6 + 100.0*ps4_6 - 200.0*ps5_6 + 100.0*ps5_7 + 100.0*ps5_5 + 
	0.25*om5_6)^2 + (-25.0*(ps5_7-ps5_5) * (ph6_6-ph4_6) + 25.0*(ps6_6-ps4_6) * 
	(ph5_7-ph5_5) - 200.0*ph5_6 + 100.0*ph6_6 + 100.0*ph4_6 - 200.0*ph5_6 + 
	100.0*ph5_7 + 100.0*ph5_5)^2 + (-200.0*om6_6 + 100.0*om7_6 + 100.0*om5_6 - 
	200.0*om6_6 + 100.0*om6_7 + 100.0*om6_5 + 1.5308084989341916e-13*ph7_6 - 
	1.5308084989341916e-13*ph5_6 - 2500.0*ph6_7 + 2500.0*ph6_5)^2 + (-200.0*ps6_6 + 
	100.0*ps7_6 + 100.0*ps5_6 - 200.0*ps6_6 + 100.0*ps6_7 + 100.0*ps6_5 + 
	0.25*om6_6)^2 + (-25.0*(ps6_7-ps6_5) * (ph7_6-ph5_6) + 25.0*(ps7_6-ps5_6) * 
	(ph6_7-ph6_5) - 200.0*ph6_6 + 100.0*ph7_6 + 100.0*ph5_6 - 200.0*ph6_6 + 
	100.0*ph6_7 + 100.0*ph6_5)^2 + (-200.0*om7_6 + 100.0*om8_6 + 100.0*om6_6 - 
	200.0*om7_6 + 100.0*om7_7 + 100.0*om7_5 + 1.5308084989341916e-13*ph8_6 - 
	1.5308084989341916e-13*ph6_6 - 2500.0*ph7_7 + 2500.0*ph7_5)^2 + (-200.0*ps7_6 + 
	100.0*ps8_6 + 100.0*ps6_6 - 200.0*ps7_6 + 100.0*ps7_7 + 100.0*ps7_5 + 
	0.25*om7_6)^2 + (-25.0*(ps7_7-ps7_5) * (ph8_6-ph6_6) + 25.0*(ps8_6-ps6_6) * 
	(ph7_7-ph7_5) - 200.0*ph7_6 + 100.0*ph8_6 + 100.0*ph6_6 - 200.0*ph7_6 + 
	100.0*ph7_7 + 100.0*ph7_5)^2 + (-200.0*om8_6 + 100.0*om9_6 + 100.0*om7_6 - 
	200.0*om8_6 + 100.0*om8_7 + 100.0*om8_5 + 1.5308084989341916e-13*ph9_6 - 
	1.5308084989341916e-13*ph7_6 - 2500.0*ph8_7 + 2500.0*ph8_5)^2 + (-200.0*ps8_6 + 
	100.0*ps9_6 + 100.0*ps7_6 - 200.0*ps8_6 + 100.0*ps8_7 + 100.0*ps8_5 + 
	0.25*om8_6)^2 + (-25.0*(ps8_7-ps8_5) * (ph9_6-ph7_6) + 25.0*(ps9_6-ps7_6) * 
	(ph8_7-ph8_5) - 200.0*ph8_6 + 100.0*ph9_6 + 100.0*ph7_6 - 200.0*ph8_6 + 
	100.0*ph8_7 + 100.0*ph8_5)^2 + (-200.0*om9_6 + 100.0*om10_6 + 100.0*om8_6 - 
	200.0*om9_6 + 100.0*om9_7 + 100.0*om9_5 + 1.5308084989341916e-13*ph10_6 - 
	1.5308084989341916e-13*ph8_6 - 2500.0*ph9_7 + 2500.0*ph9_5)^2 + (-200.0*ps9_6 + 
	100.0*ps10_6 + 100.0*ps8_6 - 200.0*ps9_6 + 100.0*ps9_7 + 100.0*ps9_5 + 
	0.25*om9_6)^2 + (-25.0*(ps9_7-ps9_5) * (ph10_6-ph8_6) + 25.0*(ps10_6-ps8_6) * 
	(ph9_7-ph9_5) - 200.0*ph9_6 + 100.0*ph10_6 + 100.0*ph8_6 - 200.0*ph9_6 + 
	100.0*ph9_7 + 100.0*ph9_5)^2 + (-200.0*omm9_7 + 100.0*omm8_7 + 100.0*omm10_7 - 
	200.0*omm9_7 + 100.0*omm9_8 + 100.0*omm9_6 + 1.5308084989341916e-13*phm8_7 - 
	1.5308084989341916e-13*phm10_7 - 2500.0*phm9_8 + 2500.0*phm9_6)^2 + 
	(-200.0*psm9_7 + 100.0*psm8_7 + 100.0*psm10_7 - 200.0*psm9_7 + 100.0*psm9_8 + 
	100.0*psm9_6 + 0.25*omm9_7)^2 + (-25.0*(psm9_8-psm9_6) * (phm8_7-phm10_7) + 
	25.0*(psm8_7-psm10_7) * (phm9_8-phm9_6) - 200.0*phm9_7 + 100.0*phm8_7 + 
	100.0*phm10_7 - 200.0*phm9_7 + 100.0*phm9_8 + 100.0*phm9_6)^2 + (-200.0*omm8_7 
	+ 100.0*omm7_7 + 100.0*omm9_7 - 200.0*omm8_7 + 100.0*omm8_8 + 100.0*omm8_6 + 
	1.5308084989341916e-13*phm7_7 - 1.5308084989341916e-13*phm9_7 - 2500.0*phm8_8 + 
	2500.0*phm8_6)^2 + (-200.0*psm8_7 + 100.0*psm7_7 + 100.0*psm9_7 - 200.0*psm8_7 
	+ 100.0*psm8_8 + 100.0*psm8_6 + 0.25*omm8_7)^2 + (-25.0*(psm8_8-psm8_6) * 
	(phm7_7-phm9_7) + 25.0*(psm7_7-psm9_7) * (phm8_8-phm8_6) - 200.0*phm8_7 + 
	100.0*phm7_7 + 100.0*phm9_7 - 200.0*phm8_7 + 100.0*phm8_8 + 100.0*phm8_6)^2 + 
	(-200.0*omm7_7 + 100.0*omm6_7 + 100.0*omm8_7 - 200.0*omm7_7 + 100.0*omm7_8 + 
	100.0*omm7_6 + 1.5308084989341916e-13*phm6_7 - 1.5308084989341916e-13*phm8_7 - 
	2500.0*phm7_8 + 2500.0*phm7_6)^2 + (-200.0*psm7_7 + 100.0*psm6_7 + 100.0*psm8_7 
	- 200.0*psm7_7 + 100.0*psm7_8 + 100.0*psm7_6 + 0.25*omm7_7)^2 + 
	(-25.0*(psm7_8-psm7_6) * (phm6_7-phm8_7) + 25.0*(psm6_7-psm8_7) * 
	(phm7_8-phm7_6) - 200.0*phm7_7 + 100.0*phm6_7 + 100.0*phm8_7 - 200.0*phm7_7 + 
	100.0*phm7_8 + 100.0*phm7_6)^2 + (-200.0*omm6_7 + 100.0*omm5_7 + 100.0*omm7_7 - 
	200.0*omm6_7 + 100.0*omm6_8 + 100.0*omm6_6 + 1.5308084989341916e-13*phm5_7 - 
	1.5308084989341916e-13*phm7_7 - 2500.0*phm6_8 + 2500.0*phm6_6)^2 + 
	(-200.0*psm6_7 + 100.0*psm5_7 + 100.0*psm7_7 - 200.0*psm6_7 + 100.0*psm6_8 + 
	100.0*psm6_6 + 0.25*omm6_7)^2 + (-25.0*(psm6_8-psm6_6) * (phm5_7-phm7_7) + 
	25.0*(psm5_7-psm7_7) * (phm6_8-phm6_6) - 200.0*phm6_7 + 100.0*phm5_7 + 
	100.0*phm7_7 - 200.0*phm6_7 + 100.0*phm6_8 + 100.0*phm6_6)^2 + (-200.0*omm5_7 + 
	100.0*omm4_7 + 100.0*omm6_7 - 200.0*omm5_7 + 100.0*omm5_8 + 100.0*omm5_6 + 
	1.5308084989341916e-13*phm4_7 - 1.5308084989341916e-13*phm6_7 - 2500.0*phm5_8 + 
	2500.0*phm5_6)^2 + (-200.0*psm5_7 + 100.0*psm4_7 + 100.0*psm6_7 - 200.0*psm5_7 
	+ 100.0*psm5_8 + 100.0*psm5_6 + 0.25*omm5_7)^2 + (-25.0*(psm5_8-psm5_6) * 
	(phm4_7-phm6_7) + 25.0*(psm4_7-psm6_7) * (phm5_8-phm5_6) - 200.0*phm5_7 + 
	100.0*phm4_7 + 100.0*phm6_7 - 200.0*phm5_7 + 100.0*phm5_8 + 100.0*phm5_6)^2 + 
	(-200.0*omm4_7 + 100.0*omm3_7 + 100.0*omm5_7 - 200.0*omm4_7 + 100.0*omm4_8 + 
	100.0*omm4_6 + 1.5308084989341916e-13*phm3_7 - 1.5308084989341916e-13*phm5_7 - 
	2500.0*phm4_8 + 2500.0*phm4_6)^2 + (-200.0*psm4_7 + 100.0*psm3_7 + 100.0*psm5_7 
	- 200.0*psm4_7 + 100.0*psm4_8 + 100.0*psm4_6 + 0.25*omm4_7)^2 + 
	(-25.0*(psm4_8-psm4_6) * (phm3_7-phm5_7) + 25.0*(psm3_7-psm5_7) * 
	(phm4_8-phm4_6) - 200.0*phm4_7 + 100.0*phm3_7 + 100.0*phm5_7 - 200.0*phm4_7 + 
	100.0*phm4_8 + 100.0*phm4_6)^2 + (-200.0*omm3_7 + 100.0*omm2_7 + 100.0*omm4_7 - 
	200.0*omm3_7 + 100.0*omm3_8 + 100.0*omm3_6 + 1.5308084989341916e-13*phm2_7 - 
	1.5308084989341916e-13*phm4_7 - 2500.0*phm3_8 + 2500.0*phm3_6)^2 + 
	(-200.0*psm3_7 + 100.0*psm2_7 + 100.0*psm4_7 - 200.0*psm3_7 + 100.0*psm3_8 + 
	100.0*psm3_6 + 0.25*omm3_7)^2 + (-25.0*(psm3_8-psm3_6) * (phm2_7-phm4_7) + 
	25.0*(psm2_7-psm4_7) * (phm3_8-phm3_6) - 200.0*phm3_7 + 100.0*phm2_7 + 
	100.0*phm4_7 - 200.0*phm3_7 + 100.0*phm3_8 + 100.0*phm3_6)^2 + (-200.0*omm2_7 + 
	100.0*omm1_7 + 100.0*omm3_7 - 200.0*omm2_7 + 100.0*omm2_8 + 100.0*omm2_6 + 
	1.5308084989341916e-13*phm1_7 - 1.5308084989341916e-13*phm3_7 - 2500.0*phm2_8 + 
	2500.0*phm2_6)^2 + (-200.0*psm2_7 + 100.0*psm1_7 + 100.0*psm3_7 - 200.0*psm2_7 
	+ 100.0*psm2_8 + 100.0*psm2_6 + 0.25*omm2_7)^2 + (-25.0*(psm2_8-psm2_6) * 
	(phm1_7-phm3_7) + 25.0*(psm1_7-psm3_7) * (phm2_8-phm2_6) - 200.0*phm2_7 + 
	100.0*phm1_7 + 100.0*phm3_7 - 200.0*phm2_7 + 100.0*phm2_8 + 100.0*phm2_6)^2 + 
	(-200.0*omm1_7 + 100.0*om0_7 + 100.0*omm2_7 - 200.0*omm1_7 + 100.0*omm1_8 + 
	100.0*omm1_6 + 1.5308084989341916e-13*ph0_7 - 1.5308084989341916e-13*phm2_7 - 
	2500.0*phm1_8 + 2500.0*phm1_6)^2 + (-200.0*psm1_7 + 100.0*ps0_7 + 100.0*psm2_7 
	- 200.0*psm1_7 + 100.0*psm1_8 + 100.0*psm1_6 + 0.25*omm1_7)^2 + 
	(-25.0*(psm1_8-psm1_6) * (ph0_7-phm2_7) + 25.0*(ps0_7-psm2_7) * (phm1_8-phm1_6) 
	- 200.0*phm1_7 + 100.0*ph0_7 + 100.0*phm2_7 - 200.0*phm1_7 + 100.0*phm1_8 + 
	100.0*phm1_6)^2 + (-200.0*om0_7 + 100.0*om1_7 + 100.0*omm1_7 - 200.0*om0_7 + 
	100.0*om0_8 + 100.0*om0_6 + 1.5308084989341916e-13*ph1_7 - 
	1.5308084989341916e-13*phm1_7 - 2500.0*ph0_8 + 2500.0*ph0_6)^2 + (-200.0*ps0_7 
	+ 100.0*ps1_7 + 100.0*psm1_7 - 200.0*ps0_7 + 100.0*ps0_8 + 100.0*ps0_6 + 
	0.25*om0_7)^2 + (-25.0*(ps0_8-ps0_6) * (ph1_7-phm1_7) + 25.0*(ps1_7-psm1_7) * 
	(ph0_8-ph0_6) - 200.0*ph0_7 + 100.0*ph1_7 + 100.0*phm1_7 - 200.0*ph0_7 + 
	100.0*ph0_8 + 100.0*ph0_6)^2 + (-200.0*om1_7 + 100.0*om2_7 + 100.0*om0_7 - 
	200.0*om1_7 + 100.0*om1_8 + 100.0*om1_6 + 1.5308084989341916e-13*ph2_7 - 
	1.5308084989341916e-13*ph0_7 - 2500.0*ph1_8 + 2500.0*ph1_6)^2 + (-200.0*ps1_7 + 
	100.0*ps2_7 + 100.0*ps0_7 - 200.0*ps1_7 + 100.0*ps1_8 + 100.0*ps1_6 + 
	0.25*om1_7)^2 + (-25.0*(ps1_8-ps1_6) * (ph2_7-ph0_7) + 25.0*(ps2_7-ps0_7) * 
	(ph1_8-ph1_6) - 200.0*ph1_7 + 100.0*ph2_7 + 100.0*ph0_7 - 200.0*ph1_7 + 
	100.0*ph1_8 + 100.0*ph1_6)^2 + (-200.0*om2_7 + 100.0*om3_7 + 100.0*om1_7 - 
	200.0*om2_7 + 100.0*om2_8 + 100.0*om2_6 + 1.5308084989341916e-13*ph3_7 - 
	1.5308084989341916e-13*ph1_7 - 2500.0*ph2_8 + 2500.0*ph2_6)^2 + (-200.0*ps2_7 + 
	100.0*ps3_7 + 100.0*ps1_7 - 200.0*ps2_7 + 100.0*ps2_8 + 100.0*ps2_6 + 
	0.25*om2_7)^2 + (-25.0*(ps2_8-ps2_6) * (ph3_7-ph1_7) + 25.0*(ps3_7-ps1_7) * 
	(ph2_8-ph2_6) - 200.0*ph2_7 + 100.0*ph3_7 + 100.0*ph1_7 - 200.0*ph2_7 + 
	100.0*ph2_8 + 100.0*ph2_6)^2 + (-200.0*om3_7 + 100.0*om4_7 + 100.0*om2_7 - 
	200.0*om3_7 + 100.0*om3_8 + 100.0*om3_6 + 1.5308084989341916e-13*ph4_7 - 
	1.5308084989341916e-13*ph2_7 - 2500.0*ph3_8 + 2500.0*ph3_6)^2 + (-200.0*ps3_7 + 
	100.0*ps4_7 + 100.0*ps2_7 - 200.0*ps3_7 + 100.0*ps3_8 + 100.0*ps3_6 + 
	0.25*om3_7)^2 + (-25.0*(ps3_8-ps3_6) * (ph4_7-ph2_7) + 25.0*(ps4_7-ps2_7) * 
	(ph3_8-ph3_6) - 200.0*ph3_7 + 100.0*ph4_7 + 100.0*ph2_7 - 200.0*ph3_7 + 
	100.0*ph3_8 + 100.0*ph3_6)^2 + (-200.0*om4_7 + 100.0*om5_7 + 100.0*om3_7 - 
	200.0*om4_7 + 100.0*om4_8 + 100.0*om4_6 + 1.5308084989341916e-13*ph5_7 - 
	1.5308084989341916e-13*ph3_7 - 2500.0*ph4_8 + 2500.0*ph4_6)^2 + (-200.0*ps4_7 + 
	100.0*ps5_7 + 100.0*ps3_7 - 200.0*ps4_7 + 100.0*ps4_8 + 100.0*ps4_6 + 
	0.25*om4_7)^2 + (-25.0*(ps4_8-ps4_6) * (ph5_7-ph3_7) + 25.0*(ps5_7-ps3_7) * 
	(ph4_8-ph4_6) - 200.0*ph4_7 + 100.0*ph5_7 + 100.0*ph3_7 - 200.0*ph4_7 + 
	100.0*ph4_8 + 100.0*ph4_6)^2 + (-200.0*om5_7 + 100.0*om6_7 + 100.0*om4_7 - 
	200.0*om5_7 + 100.0*om5_8 + 100.0*om5_6 + 1.5308084989341916e-13*ph6_7 - 
	1.5308084989341916e-13*ph4_7 - 2500.0*ph5_8 + 2500.0*ph5_6)^2 + (-200.0*ps5_7 + 
	100.0*ps6_7 + 100.0*ps4_7 - 200.0*ps5_7 + 100.0*ps5_8 + 100.0*ps5_6 + 
	0.25*om5_7)^2 + (-25.0*(ps5_8-ps5_6) * (ph6_7-ph4_7) + 25.0*(ps6_7-ps4_7) * 
	(ph5_8-ph5_6) - 200.0*ph5_7 + 100.0*ph6_7 + 100.0*ph4_7 - 200.0*ph5_7 + 
	100.0*ph5_8 + 100.0*ph5_6)^2 + (-200.0*om6_7 + 100.0*om7_7 + 100.0*om5_7 - 
	200.0*om6_7 + 100.0*om6_8 + 100.0*om6_6 + 1.5308084989341916e-13*ph7_7 - 
	1.5308084989341916e-13*ph5_7 - 2500.0*ph6_8 + 2500.0*ph6_6)^2 + (-200.0*ps6_7 + 
	100.0*ps7_7 + 100.0*ps5_7 - 200.0*ps6_7 + 100.0*ps6_8 + 100.0*ps6_6 + 
	0.25*om6_7)^2 + (-25.0*(ps6_8-ps6_6) * (ph7_7-ph5_7) + 25.0*(ps7_7-ps5_7) * 
	(ph6_8-ph6_6) - 200.0*ph6_7 + 100.0*ph7_7 + 100.0*ph5_7 - 200.0*ph6_7 + 
	100.0*ph6_8 + 100.0*ph6_6)^2 + (-200.0*om7_7 + 100.0*om8_7 + 100.0*om6_7 - 
	200.0*om7_7 + 100.0*om7_8 + 100.0*om7_6 + 1.5308084989341916e-13*ph8_7 - 
	1.5308084989341916e-13*ph6_7 - 2500.0*ph7_8 + 2500.0*ph7_6)^2 + (-200.0*ps7_7 + 
	100.0*ps8_7 + 100.0*ps6_7 - 200.0*ps7_7 + 100.0*ps7_8 + 100.0*ps7_6 + 
	0.25*om7_7)^2 + (-25.0*(ps7_8-ps7_6) * (ph8_7-ph6_7) + 25.0*(ps8_7-ps6_7) * 
	(ph7_8-ph7_6) - 200.0*ph7_7 + 100.0*ph8_7 + 100.0*ph6_7 - 200.0*ph7_7 + 
	100.0*ph7_8 + 100.0*ph7_6)^2 + (-200.0*om8_7 + 100.0*om9_7 + 100.0*om7_7 - 
	200.0*om8_7 + 100.0*om8_8 + 100.0*om8_6 + 1.5308084989341916e-13*ph9_7 - 
	1.5308084989341916e-13*ph7_7 - 2500.0*ph8_8 + 2500.0*ph8_6)^2 + (-200.0*ps8_7 + 
	100.0*ps9_7 + 100.0*ps7_7 - 200.0*ps8_7 + 100.0*ps8_8 + 100.0*ps8_6 + 
	0.25*om8_7)^2 + (-25.0*(ps8_8-ps8_6) * (ph9_7-ph7_7) + 25.0*(ps9_7-ps7_7) * 
	(ph8_8-ph8_6) - 200.0*ph8_7 + 100.0*ph9_7 + 100.0*ph7_7 - 200.0*ph8_7 + 
	100.0*ph8_8 + 100.0*ph8_6)^2 + (-200.0*om9_7 + 100.0*om10_7 + 100.0*om8_7 - 
	200.0*om9_7 + 100.0*om9_8 + 100.0*om9_6 + 1.5308084989341916e-13*ph10_7 - 
	1.5308084989341916e-13*ph8_7 - 2500.0*ph9_8 + 2500.0*ph9_6)^2 + (-200.0*ps9_7 + 
	100.0*ps10_7 + 100.0*ps8_7 - 200.0*ps9_7 + 100.0*ps9_8 + 100.0*ps9_6 + 
	0.25*om9_7)^2 + (-25.0*(ps9_8-ps9_6) * (ph10_7-ph8_7) + 25.0*(ps10_7-ps8_7) * 
	(ph9_8-ph9_6) - 200.0*ph9_7 + 100.0*ph10_7 + 100.0*ph8_7 - 200.0*ph9_7 + 
	100.0*ph9_8 + 100.0*ph9_6)^2 + (-200.0*omm9_8 + 100.0*omm8_8 + 100.0*omm10_8 - 
	200.0*omm9_8 + 100.0*omm9_9 + 100.0*omm9_7 + 1.5308084989341916e-13*phm8_8 - 
	1.5308084989341916e-13*phm10_8 - 2500.0*phm9_9 + 2500.0*phm9_7)^2 + 
	(-200.0*psm9_8 + 100.0*psm8_8 + 100.0*psm10_8 - 200.0*psm9_8 + 100.0*psm9_9 + 
	100.0*psm9_7 + 0.25*omm9_8)^2 + (-25.0*(psm9_9-psm9_7) * (phm8_8-phm10_8) + 
	25.0*(psm8_8-psm10_8) * (phm9_9-phm9_7) - 200.0*phm9_8 + 100.0*phm8_8 + 
	100.0*phm10_8 - 200.0*phm9_8 + 100.0*phm9_9 + 100.0*phm9_7)^2 + (-200.0*omm8_8 
	+ 100.0*omm7_8 + 100.0*omm9_8 - 200.0*omm8_8 + 100.0*omm8_9 + 100.0*omm8_7 + 
	1.5308084989341916e-13*phm7_8 - 1.5308084989341916e-13*phm9_8 - 2500.0*phm8_9 + 
	2500.0*phm8_7)^2 + (-200.0*psm8_8 + 100.0*psm7_8 + 100.0*psm9_8 - 200.0*psm8_8 
	+ 100.0*psm8_9 + 100.0*psm8_7 + 0.25*omm8_8)^2 + (-25.0*(psm8_9-psm8_7) * 
	(phm7_8-phm9_8) + 25.0*(psm7_8-psm9_8) * (phm8_9-phm8_7) - 200.0*phm8_8 + 
	100.0*phm7_8 + 100.0*phm9_8 - 200.0*phm8_8 + 100.0*phm8_9 + 100.0*phm8_7)^2 + 
	(-200.0*omm7_8 + 100.0*omm6_8 + 100.0*omm8_8 - 200.0*omm7_8 + 100.0*omm7_9 + 
	100.0*omm7_7 + 1.5308084989341916e-13*phm6_8 - 1.5308084989341916e-13*phm8_8 - 
	2500.0*phm7_9 + 2500.0*phm7_7)^2 + (-200.0*psm7_8 + 100.0*psm6_8 + 100.0*psm8_8 
	- 200.0*psm7_8 + 100.0*psm7_9 + 100.0*psm7_7 + 0.25*omm7_8)^2 + 
	(-25.0*(psm7_9-psm7_7) * (phm6_8-phm8_8) + 25.0*(psm6_8-psm8_8) * 
	(phm7_9-phm7_7) - 200.0*phm7_8 + 100.0*phm6_8 + 100.0*phm8_8 - 200.0*phm7_8 + 
	100.0*phm7_9 + 100.0*phm7_7)^2 + (-200.0*omm6_8 + 100.0*omm5_8 + 100.0*omm7_8 - 
	200.0*omm6_8 + 100.0*omm6_9 + 100.0*omm6_7 + 1.5308084989341916e-13*phm5_8 - 
	1.5308084989341916e-13*phm7_8 - 2500.0*phm6_9 + 2500.0*phm6_7)^2 + 
	(-200.0*psm6_8 + 100.0*psm5_8 + 100.0*psm7_8 - 200.0*psm6_8 + 100.0*psm6_9 + 
	100.0*psm6_7 + 0.25*omm6_8)^2 + (-25.0*(psm6_9-psm6_7) * (phm5_8-phm7_8) + 
	25.0*(psm5_8-psm7_8) * (phm6_9-phm6_7) - 200.0*phm6_8 + 100.0*phm5_8 + 
	100.0*phm7_8 - 200.0*phm6_8 + 100.0*phm6_9 + 100.0*phm6_7)^2 + (-200.0*omm5_8 + 
	100.0*omm4_8 + 100.0*omm6_8 - 200.0*omm5_8 + 100.0*omm5_9 + 100.0*omm5_7 + 
	1.5308084989341916e-13*phm4_8 - 1.5308084989341916e-13*phm6_8 - 2500.0*phm5_9 + 
	2500.0*phm5_7)^2 + (-200.0*psm5_8 + 100.0*psm4_8 + 100.0*psm6_8 - 200.0*psm5_8 
	+ 100.0*psm5_9 + 100.0*psm5_7 + 0.25*omm5_8)^2 + (-25.0*(psm5_9-psm5_7) * 
	(phm4_8-phm6_8) + 25.0*(psm4_8-psm6_8) * (phm5_9-phm5_7) - 200.0*phm5_8 + 
	100.0*phm4_8 + 100.0*phm6_8 - 200.0*phm5_8 + 100.0*phm5_9 + 100.0*phm5_7)^2 + 
	(-200.0*omm4_8 + 100.0*omm3_8 + 100.0*omm5_8 - 200.0*omm4_8 + 100.0*omm4_9 + 
	100.0*omm4_7 + 1.5308084989341916e-13*phm3_8 - 1.5308084989341916e-13*phm5_8 - 
	2500.0*phm4_9 + 2500.0*phm4_7)^2 + (-200.0*psm4_8 + 100.0*psm3_8 + 100.0*psm5_8 
	- 200.0*psm4_8 + 100.0*psm4_9 + 100.0*psm4_7 + 0.25*omm4_8)^2 + 
	(-25.0*(psm4_9-psm4_7) * (phm3_8-phm5_8) + 25.0*(psm3_8-psm5_8) * 
	(phm4_9-phm4_7) - 200.0*phm4_8 + 100.0*phm3_8 + 100.0*phm5_8 - 200.0*phm4_8 + 
	100.0*phm4_9 + 100.0*phm4_7)^2 + (-200.0*omm3_8 + 100.0*omm2_8 + 100.0*omm4_8 - 
	200.0*omm3_8 + 100.0*omm3_9 + 100.0*omm3_7 + 1.5308084989341916e-13*phm2_8 - 
	1.5308084989341916e-13*phm4_8 - 2500.0*phm3_9 + 2500.0*phm3_7)^2 + 
	(-200.0*psm3_8 + 100.0*psm2_8 + 100.0*psm4_8 - 200.0*psm3_8 + 100.0*psm3_9 + 
	100.0*psm3_7 + 0.25*omm3_8)^2 + (-25.0*(psm3_9-psm3_7) * (phm2_8-phm4_8) + 
	25.0*(psm2_8-psm4_8) * (phm3_9-phm3_7) - 200.0*phm3_8 + 100.0*phm2_8 + 
	100.0*phm4_8 - 200.0*phm3_8 + 100.0*phm3_9 + 100.0*phm3_7)^2 + (-200.0*omm2_8 + 
	100.0*omm1_8 + 100.0*omm3_8 - 200.0*omm2_8 + 100.0*omm2_9 + 100.0*omm2_7 + 
	1.5308084989341916e-13*phm1_8 - 1.5308084989341916e-13*phm3_8 - 2500.0*phm2_9 + 
	2500.0*phm2_7)^2 + (-200.0*psm2_8 + 100.0*psm1_8 + 100.0*psm3_8 - 200.0*psm2_8 
	+ 100.0*psm2_9 + 100.0*psm2_7 + 0.25*omm2_8)^2 + (-25.0*(psm2_9-psm2_7) * 
	(phm1_8-phm3_8) + 25.0*(psm1_8-psm3_8) * (phm2_9-phm2_7) - 200.0*phm2_8 + 
	100.0*phm1_8 + 100.0*phm3_8 - 200.0*phm2_8 + 100.0*phm2_9 + 100.0*phm2_7)^2 + 
	(-200.0*omm1_8 + 100.0*om0_8 + 100.0*omm2_8 - 200.0*omm1_8 + 100.0*omm1_9 + 
	100.0*omm1_7 + 1.5308084989341916e-13*ph0_8 - 1.5308084989341916e-13*phm2_8 - 
	2500.0*phm1_9 + 2500.0*phm1_7)^2 + (-200.0*psm1_8 + 100.0*ps0_8 + 100.0*psm2_8 
	- 200.0*psm1_8 + 100.0*psm1_9 + 100.0*psm1_7 + 0.25*omm1_8)^2 + 
	(-25.0*(psm1_9-psm1_7) * (ph0_8-phm2_8) + 25.0*(ps0_8-psm2_8) * (phm1_9-phm1_7) 
	- 200.0*phm1_8 + 100.0*ph0_8 + 100.0*phm2_8 - 200.0*phm1_8 + 100.0*phm1_9 + 
	100.0*phm1_7)^2 + (-200.0*om0_8 + 100.0*om1_8 + 100.0*omm1_8 - 200.0*om0_8 + 
	100.0*om0_9 + 100.0*om0_7 + 1.5308084989341916e-13*ph1_8 - 
	1.5308084989341916e-13*phm1_8 - 2500.0*ph0_9 + 2500.0*ph0_7)^2 + (-200.0*ps0_8 
	+ 100.0*ps1_8 + 100.0*psm1_8 - 200.0*ps0_8 + 100.0*ps0_9 + 100.0*ps0_7 + 
	0.25*om0_8)^2 + (-25.0*(ps0_9-ps0_7) * (ph1_8-phm1_8) + 25.0*(ps1_8-psm1_8) * 
	(ph0_9-ph0_7) - 200.0*ph0_8 + 100.0*ph1_8 + 100.0*phm1_8 - 200.0*ph0_8 + 
	100.0*ph0_9 + 100.0*ph0_7)^2 + (-200.0*om1_8 + 100.0*om2_8 + 100.0*om0_8 - 
	200.0*om1_8 + 100.0*om1_9 + 100.0*om1_7 + 1.5308084989341916e-13*ph2_8 - 
	1.5308084989341916e-13*ph0_8 - 2500.0*ph1_9 + 2500.0*ph1_7)^2 + (-200.0*ps1_8 + 
	100.0*ps2_8 + 100.0*ps0_8 - 200.0*ps1_8 + 100.0*ps1_9 + 100.0*ps1_7 + 
	0.25*om1_8)^2 + (-25.0*(ps1_9-ps1_7) * (ph2_8-ph0_8) + 25.0*(ps2_8-ps0_8) * 
	(ph1_9-ph1_7) - 200.0*ph1_8 + 100.0*ph2_8 + 100.0*ph0_8 - 200.0*ph1_8 + 
	100.0*ph1_9 + 100.0*ph1_7)^2 + (-200.0*om2_8 + 100.0*om3_8 + 100.0*om1_8 - 
	200.0*om2_8 + 100.0*om2_9 + 100.0*om2_7 + 1.5308084989341916e-13*ph3_8 - 
	1.5308084989341916e-13*ph1_8 - 2500.0*ph2_9 + 2500.0*ph2_7)^2 + (-200.0*ps2_8 + 
	100.0*ps3_8 + 100.0*ps1_8 - 200.0*ps2_8 + 100.0*ps2_9 + 100.0*ps2_7 + 
	0.25*om2_8)^2 + (-25.0*(ps2_9-ps2_7) * (ph3_8-ph1_8) + 25.0*(ps3_8-ps1_8) * 
	(ph2_9-ph2_7) - 200.0*ph2_8 + 100.0*ph3_8 + 100.0*ph1_8 - 200.0*ph2_8 + 
	100.0*ph2_9 + 100.0*ph2_7)^2 + (-200.0*om3_8 + 100.0*om4_8 + 100.0*om2_8 - 
	200.0*om3_8 + 100.0*om3_9 + 100.0*om3_7 + 1.5308084989341916e-13*ph4_8 - 
	1.5308084989341916e-13*ph2_8 - 2500.0*ph3_9 + 2500.0*ph3_7)^2 + (-200.0*ps3_8 + 
	100.0*ps4_8 + 100.0*ps2_8 - 200.0*ps3_8 + 100.0*ps3_9 + 100.0*ps3_7 + 
	0.25*om3_8)^2 + (-25.0*(ps3_9-ps3_7) * (ph4_8-ph2_8) + 25.0*(ps4_8-ps2_8) * 
	(ph3_9-ph3_7) - 200.0*ph3_8 + 100.0*ph4_8 + 100.0*ph2_8 - 200.0*ph3_8 + 
	100.0*ph3_9 + 100.0*ph3_7)^2 + (-200.0*om4_8 + 100.0*om5_8 + 100.0*om3_8 - 
	200.0*om4_8 + 100.0*om4_9 + 100.0*om4_7 + 1.5308084989341916e-13*ph5_8 - 
	1.5308084989341916e-13*ph3_8 - 2500.0*ph4_9 + 2500.0*ph4_7)^2 + (-200.0*ps4_8 + 
	100.0*ps5_8 + 100.0*ps3_8 - 200.0*ps4_8 + 100.0*ps4_9 + 100.0*ps4_7 + 
	0.25*om4_8)^2 + (-25.0*(ps4_9-ps4_7) * (ph5_8-ph3_8) + 25.0*(ps5_8-ps3_8) * 
	(ph4_9-ph4_7) - 200.0*ph4_8 + 100.0*ph5_8 + 100.0*ph3_8 - 200.0*ph4_8 + 
	100.0*ph4_9 + 100.0*ph4_7)^2 + (-200.0*om5_8 + 100.0*om6_8 + 100.0*om4_8 - 
	200.0*om5_8 + 100.0*om5_9 + 100.0*om5_7 + 1.5308084989341916e-13*ph6_8 - 
	1.5308084989341916e-13*ph4_8 - 2500.0*ph5_9 + 2500.0*ph5_7)^2 + (-200.0*ps5_8 + 
	100.0*ps6_8 + 100.0*ps4_8 - 200.0*ps5_8 + 100.0*ps5_9 + 100.0*ps5_7 + 
	0.25*om5_8)^2 + (-25.0*(ps5_9-ps5_7) * (ph6_8-ph4_8) + 25.0*(ps6_8-ps4_8) * 
	(ph5_9-ph5_7) - 200.0*ph5_8 + 100.0*ph6_8 + 100.0*ph4_8 - 200.0*ph5_8 + 
	100.0*ph5_9 + 100.0*ph5_7)^2 + (-200.0*om6_8 + 100.0*om7_8 + 100.0*om5_8 - 
	200.0*om6_8 + 100.0*om6_9 + 100.0*om6_7 + 1.5308084989341916e-13*ph7_8 - 
	1.5308084989341916e-13*ph5_8 - 2500.0*ph6_9 + 2500.0*ph6_7)^2 + (-200.0*ps6_8 + 
	100.0*ps7_8 + 100.0*ps5_8 - 200.0*ps6_8 + 100.0*ps6_9 + 100.0*ps6_7 + 
	0.25*om6_8)^2 + (-25.0*(ps6_9-ps6_7) * (ph7_8-ph5_8) + 25.0*(ps7_8-ps5_8) * 
	(ph6_9-ph6_7) - 200.0*ph6_8 + 100.0*ph7_8 + 100.0*ph5_8 - 200.0*ph6_8 + 
	100.0*ph6_9 + 100.0*ph6_7)^2 + (-200.0*om7_8 + 100.0*om8_8 + 100.0*om6_8 - 
	200.0*om7_8 + 100.0*om7_9 + 100.0*om7_7 + 1.5308084989341916e-13*ph8_8 - 
	1.5308084989341916e-13*ph6_8 - 2500.0*ph7_9 + 2500.0*ph7_7)^2 + (-200.0*ps7_8 + 
	100.0*ps8_8 + 100.0*ps6_8 - 200.0*ps7_8 + 100.0*ps7_9 + 100.0*ps7_7 + 
	0.25*om7_8)^2 + (-25.0*(ps7_9-ps7_7) * (ph8_8-ph6_8) + 25.0*(ps8_8-ps6_8) * 
	(ph7_9-ph7_7) - 200.0*ph7_8 + 100.0*ph8_8 + 100.0*ph6_8 - 200.0*ph7_8 + 
	100.0*ph7_9 + 100.0*ph7_7)^2 + (-200.0*om8_8 + 100.0*om9_8 + 100.0*om7_8 - 
	200.0*om8_8 + 100.0*om8_9 + 100.0*om8_7 + 1.5308084989341916e-13*ph9_8 - 
	1.5308084989341916e-13*ph7_8 - 2500.0*ph8_9 + 2500.0*ph8_7)^2 + (-200.0*ps8_8 + 
	100.0*ps9_8 + 100.0*ps7_8 - 200.0*ps8_8 + 100.0*ps8_9 + 100.0*ps8_7 + 
	0.25*om8_8)^2 + (-25.0*(ps8_9-ps8_7) * (ph9_8-ph7_8) + 25.0*(ps9_8-ps7_8) * 
	(ph8_9-ph8_7) - 200.0*ph8_8 + 100.0*ph9_8 + 100.0*ph7_8 - 200.0*ph8_8 + 
	100.0*ph8_9 + 100.0*ph8_7)^2 + (-200.0*om9_8 + 100.0*om10_8 + 100.0*om8_8 - 
	200.0*om9_8 + 100.0*om9_9 + 100.0*om9_7 + 1.5308084989341916e-13*ph10_8 - 
	1.5308084989341916e-13*ph8_8 - 2500.0*ph9_9 + 2500.0*ph9_7)^2 + (-200.0*ps9_8 + 
	100.0*ps10_8 + 100.0*ps8_8 - 200.0*ps9_8 + 100.0*ps9_9 + 100.0*ps9_7 + 
	0.25*om9_8)^2 + (-25.0*(ps9_9-ps9_7) * (ph10_8-ph8_8) + 25.0*(ps10_8-ps8_8) * 
	(ph9_9-ph9_7) - 200.0*ph9_8 + 100.0*ph10_8 + 100.0*ph8_8 - 200.0*ph9_8 + 
	100.0*ph9_9 + 100.0*ph9_7)^2 + (-200.0*omm9_9 + 100.0*omm8_9 + 100.0*omm10_9 - 
	200.0*omm9_9 + 100.0*omm9_10 + 100.0*omm9_8 + 1.5308084989341916e-13*phm8_9 - 
	1.5308084989341916e-13*phm10_9 - 2500.0*phm9_10 + 2500.0*phm9_8)^2 + 
	(-200.0*psm9_9 + 100.0*psm8_9 + 100.0*psm10_9 - 200.0*psm9_9 + 100.0*psm9_10 + 
	100.0*psm9_8 + 0.25*omm9_9)^2 + (-25.0*(psm9_10-psm9_8) * (phm8_9-phm10_9) + 
	25.0*(psm8_9-psm10_9) * (phm9_10-phm9_8) - 200.0*phm9_9 + 100.0*phm8_9 + 
	100.0*phm10_9 - 200.0*phm9_9 + 100.0*phm9_10 + 100.0*phm9_8)^2 + (-200.0*omm8_9 
	+ 100.0*omm7_9 + 100.0*omm9_9 - 200.0*omm8_9 + 100.0*omm8_10 + 100.0*omm8_8 + 
	1.5308084989341916e-13*phm7_9 - 1.5308084989341916e-13*phm9_9 - 2500.0*phm8_10 
	+ 2500.0*phm8_8)^2 + (-200.0*psm8_9 + 100.0*psm7_9 + 100.0*psm9_9 - 
	200.0*psm8_9 + 100.0*psm8_10 + 100.0*psm8_8 + 0.25*omm8_9)^2 + 
	(-25.0*(psm8_10-psm8_8) * (phm7_9-phm9_9) + 25.0*(psm7_9-psm9_9) * 
	(phm8_10-phm8_8) - 200.0*phm8_9 + 100.0*phm7_9 + 100.0*phm9_9 - 200.0*phm8_9 + 
	100.0*phm8_10 + 100.0*phm8_8)^2 + (-200.0*omm7_9 + 100.0*omm6_9 + 100.0*omm8_9 
	- 200.0*omm7_9 + 100.0*omm7_10 + 100.0*omm7_8 + 1.5308084989341916e-13*phm6_9 - 
	1.5308084989341916e-13*phm8_9 - 2500.0*phm7_10 + 2500.0*phm7_8)^2 + 
	(-200.0*psm7_9 + 100.0*psm6_9 + 100.0*psm8_9 - 200.0*psm7_9 + 100.0*psm7_10 + 
	100.0*psm7_8 + 0.25*omm7_9)^2 + (-25.0*(psm7_10-psm7_8) * (phm6_9-phm8_9) + 
	25.0*(psm6_9-psm8_9) * (phm7_10-phm7_8) - 200.0*phm7_9 + 100.0*phm6_9 + 
	100.0*phm8_9 - 200.0*phm7_9 + 100.0*phm7_10 + 100.0*phm7_8)^2 + (-200.0*omm6_9 
	+ 100.0*omm5_9 + 100.0*omm7_9 - 200.0*omm6_9 + 100.0*omm6_10 + 100.0*omm6_8 + 
	1.5308084989341916e-13*phm5_9 - 1.5308084989341916e-13*phm7_9 - 2500.0*phm6_10 
	+ 2500.0*phm6_8)^2 + (-200.0*psm6_9 + 100.0*psm5_9 + 100.0*psm7_9 - 
	200.0*psm6_9 + 100.0*psm6_10 + 100.0*psm6_8 + 0.25*omm6_9)^2 + 
	(-25.0*(psm6_10-psm6_8) * (phm5_9-phm7_9) + 25.0*(psm5_9-psm7_9) * 
	(phm6_10-phm6_8) - 200.0*phm6_9 + 100.0*phm5_9 + 100.0*phm7_9 - 200.0*phm6_9 + 
	100.0*phm6_10 + 100.0*phm6_8)^2 + (-200.0*omm5_9 + 100.0*omm4_9 + 100.0*omm6_9 
	- 200.0*omm5_9 + 100.0*omm5_10 + 100.0*omm5_8 + 1.5308084989341916e-13*phm4_9 - 
	1.5308084989341916e-13*phm6_9 - 2500.0*phm5_10 + 2500.0*phm5_8)^2 + 
	(-200.0*psm5_9 + 100.0*psm4_9 + 100.0*psm6_9 - 200.0*psm5_9 + 100.0*psm5_10 + 
	100.0*psm5_8 + 0.25*omm5_9)^2 + (-25.0*(psm5_10-psm5_8) * (phm4_9-phm6_9) + 
	25.0*(psm4_9-psm6_9) * (phm5_10-phm5_8) - 200.0*phm5_9 + 100.0*phm4_9 + 
	100.0*phm6_9 - 200.0*phm5_9 + 100.0*phm5_10 + 100.0*phm5_8)^2 + (-200.0*omm4_9 
	+ 100.0*omm3_9 + 100.0*omm5_9 - 200.0*omm4_9 + 100.0*omm4_10 + 100.0*omm4_8 + 
	1.5308084989341916e-13*phm3_9 - 1.5308084989341916e-13*phm5_9 - 2500.0*phm4_10 
	+ 2500.0*phm4_8)^2 + (-200.0*psm4_9 + 100.0*psm3_9 + 100.0*psm5_9 - 
	200.0*psm4_9 + 100.0*psm4_10 + 100.0*psm4_8 + 0.25*omm4_9)^2 + 
	(-25.0*(psm4_10-psm4_8) * (phm3_9-phm5_9) + 25.0*(psm3_9-psm5_9) * 
	(phm4_10-phm4_8) - 200.0*phm4_9 + 100.0*phm3_9 + 100.0*phm5_9 - 200.0*phm4_9 + 
	100.0*phm4_10 + 100.0*phm4_8)^2 + (-200.0*omm3_9 + 100.0*omm2_9 + 100.0*omm4_9 
	- 200.0*omm3_9 + 100.0*omm3_10 + 100.0*omm3_8 + 1.5308084989341916e-13*phm2_9 - 
	1.5308084989341916e-13*phm4_9 - 2500.0*phm3_10 + 2500.0*phm3_8)^2 + 
	(-200.0*psm3_9 + 100.0*psm2_9 + 100.0*psm4_9 - 200.0*psm3_9 + 100.0*psm3_10 + 
	100.0*psm3_8 + 0.25*omm3_9)^2 + (-25.0*(psm3_10-psm3_8) * (phm2_9-phm4_9) + 
	25.0*(psm2_9-psm4_9) * (phm3_10-phm3_8) - 200.0*phm3_9 + 100.0*phm2_9 + 
	100.0*phm4_9 - 200.0*phm3_9 + 100.0*phm3_10 + 100.0*phm3_8)^2 + (-200.0*omm2_9 
	+ 100.0*omm1_9 + 100.0*omm3_9 - 200.0*omm2_9 + 100.0*omm2_10 + 100.0*omm2_8 + 
	1.5308084989341916e-13*phm1_9 - 1.5308084989341916e-13*phm3_9 - 2500.0*phm2_10 
	+ 2500.0*phm2_8)^2 + (-200.0*psm2_9 + 100.0*psm1_9 + 100.0*psm3_9 - 
	200.0*psm2_9 + 100.0*psm2_10 + 100.0*psm2_8 + 0.25*omm2_9)^2 + 
	(-25.0*(psm2_10-psm2_8) * (phm1_9-phm3_9) + 25.0*(psm1_9-psm3_9) * 
	(phm2_10-phm2_8) - 200.0*phm2_9 + 100.0*phm1_9 + 100.0*phm3_9 - 200.0*phm2_9 + 
	100.0*phm2_10 + 100.0*phm2_8)^2 + (-200.0*omm1_9 + 100.0*om0_9 + 100.0*omm2_9 - 
	200.0*omm1_9 + 100.0*omm1_10 + 100.0*omm1_8 + 1.5308084989341916e-13*ph0_9 - 
	1.5308084989341916e-13*phm2_9 - 2500.0*phm1_10 + 2500.0*phm1_8)^2 + 
	(-200.0*psm1_9 + 100.0*ps0_9 + 100.0*psm2_9 - 200.0*psm1_9 + 100.0*psm1_10 + 
	100.0*psm1_8 + 0.25*omm1_9)^2 + (-25.0*(psm1_10-psm1_8) * (ph0_9-phm2_9) + 
	25.0*(ps0_9-psm2_9) * (phm1_10-phm1_8) - 200.0*phm1_9 + 100.0*ph0_9 + 
	100.0*phm2_9 - 200.0*phm1_9 + 100.0*phm1_10 + 100.0*phm1_8)^2 + (-200.0*om0_9 + 
	100.0*om1_9 + 100.0*omm1_9 - 200.0*om0_9 + 100.0*om0_10 + 100.0*om0_8 + 
	1.5308084989341916e-13*ph1_9 - 1.5308084989341916e-13*phm1_9 - 2500.0*ph0_10 + 
	2500.0*ph0_8)^2 + (-200.0*ps0_9 + 100.0*ps1_9 + 100.0*psm1_9 - 200.0*ps0_9 + 
	100.0*ps0_10 + 100.0*ps0_8 + 0.25*om0_9)^2 + (-25.0*(ps0_10-ps0_8) * 
	(ph1_9-phm1_9) + 25.0*(ps1_9-psm1_9) * (ph0_10-ph0_8) - 200.0*ph0_9 + 
	100.0*ph1_9 + 100.0*phm1_9 - 200.0*ph0_9 + 100.0*ph0_10 + 100.0*ph0_8)^2 + 
	(-200.0*om1_9 + 100.0*om2_9 + 100.0*om0_9 - 200.0*om1_9 + 100.0*om1_10 + 
	100.0*om1_8 + 1.5308084989341916e-13*ph2_9 - 1.5308084989341916e-13*ph0_9 - 
	2500.0*ph1_10 + 2500.0*ph1_8)^2 + (-200.0*ps1_9 + 100.0*ps2_9 + 100.0*ps0_9 - 
	200.0*ps1_9 + 100.0*ps1_10 + 100.0*ps1_8 + 0.25*om1_9)^2 + 
	(-25.0*(ps1_10-ps1_8) * (ph2_9-ph0_9) + 25.0*(ps2_9-ps0_9) * (ph1_10-ph1_8) - 
	200.0*ph1_9 + 100.0*ph2_9 + 100.0*ph0_9 - 200.0*ph1_9 + 100.0*ph1_10 + 
	100.0*ph1_8)^2 + (-200.0*om2_9 + 100.0*om3_9 + 100.0*om1_9 - 200.0*om2_9 + 
	100.0*om2_10 + 100.0*om2_8 + 1.5308084989341916e-13*ph3_9 - 
	1.5308084989341916e-13*ph1_9 - 2500.0*ph2_10 + 2500.0*ph2_8)^2 + (-200.0*ps2_9 
	+ 100.0*ps3_9 + 100.0*ps1_9 - 200.0*ps2_9 + 100.0*ps2_10 + 100.0*ps2_8 + 
	0.25*om2_9)^2 + (-25.0*(ps2_10-ps2_8) * (ph3_9-ph1_9) + 25.0*(ps3_9-ps1_9) * 
	(ph2_10-ph2_8) - 200.0*ph2_9 + 100.0*ph3_9 + 100.0*ph1_9 - 200.0*ph2_9 + 
	100.0*ph2_10 + 100.0*ph2_8)^2 + (-200.0*om3_9 + 100.0*om4_9 + 100.0*om2_9 - 
	200.0*om3_9 + 100.0*om3_10 + 100.0*om3_8 + 1.5308084989341916e-13*ph4_9 - 
	1.5308084989341916e-13*ph2_9 - 2500.0*ph3_10 + 2500.0*ph3_8)^2 + (-200.0*ps3_9 
	+ 100.0*ps4_9 + 100.0*ps2_9 - 200.0*ps3_9 + 100.0*ps3_10 + 100.0*ps3_8 + 
	0.25*om3_9)^2 + (-25.0*(ps3_10-ps3_8) * (ph4_9-ph2_9) + 25.0*(ps4_9-ps2_9) * 
	(ph3_10-ph3_8) - 200.0*ph3_9 + 100.0*ph4_9 + 100.0*ph2_9 - 200.0*ph3_9 + 
	100.0*ph3_10 + 100.0*ph3_8)^2 + (-200.0*om4_9 + 100.0*om5_9 + 100.0*om3_9 - 
	200.0*om4_9 + 100.0*om4_10 + 100.0*om4_8 + 1.5308084989341916e-13*ph5_9 - 
	1.5308084989341916e-13*ph3_9 - 2500.0*ph4_10 + 2500.0*ph4_8)^2 + (-200.0*ps4_9 
	+ 100.0*ps5_9 + 100.0*ps3_9 - 200.0*ps4_9 + 100.0*ps4_10 + 100.0*ps4_8 + 
	0.25*om4_9)^2 + (-25.0*(ps4_10-ps4_8) * (ph5_9-ph3_9) + 25.0*(ps5_9-ps3_9) * 
	(ph4_10-ph4_8) - 200.0*ph4_9 + 100.0*ph5_9 + 100.0*ph3_9 - 200.0*ph4_9 + 
	100.0*ph4_10 + 100.0*ph4_8)^2 + (-200.0*om5_9 + 100.0*om6_9 + 100.0*om4_9 - 
	200.0*om5_9 + 100.0*om5_10 + 100.0*om5_8 + 1.5308084989341916e-13*ph6_9 - 
	1.5308084989341916e-13*ph4_9 - 2500.0*ph5_10 + 2500.0*ph5_8)^2 + (-200.0*ps5_9 
	+ 100.0*ps6_9 + 100.0*ps4_9 - 200.0*ps5_9 + 100.0*ps5_10 + 100.0*ps5_8 + 
	0.25*om5_9)^2 + (-25.0*(ps5_10-ps5_8) * (ph6_9-ph4_9) + 25.0*(ps6_9-ps4_9) * 
	(ph5_10-ph5_8) - 200.0*ph5_9 + 100.0*ph6_9 + 100.0*ph4_9 - 200.0*ph5_9 + 
	100.0*ph5_10 + 100.0*ph5_8)^2 + (-200.0*om6_9 + 100.0*om7_9 + 100.0*om5_9 - 
	200.0*om6_9 + 100.0*om6_10 + 100.0*om6_8 + 1.5308084989341916e-13*ph7_9 - 
	1.5308084989341916e-13*ph5_9 - 2500.0*ph6_10 + 2500.0*ph6_8)^2 + (-200.0*ps6_9 
	+ 100.0*ps7_9 + 100.0*ps5_9 - 200.0*ps6_9 + 100.0*ps6_10 + 100.0*ps6_8 + 
	0.25*om6_9)^2 + (-25.0*(ps6_10-ps6_8) * (ph7_9-ph5_9) + 25.0*(ps7_9-ps5_9) * 
	(ph6_10-ph6_8) - 200.0*ph6_9 + 100.0*ph7_9 + 100.0*ph5_9 - 200.0*ph6_9 + 
	100.0*ph6_10 + 100.0*ph6_8)^2 + (-200.0*om7_9 + 100.0*om8_9 + 100.0*om6_9 - 
	200.0*om7_9 + 100.0*om7_10 + 100.0*om7_8 + 1.5308084989341916e-13*ph8_9 - 
	1.5308084989341916e-13*ph6_9 - 2500.0*ph7_10 + 2500.0*ph7_8)^2 + (-200.0*ps7_9 
	+ 100.0*ps8_9 + 100.0*ps6_9 - 200.0*ps7_9 + 100.0*ps7_10 + 100.0*ps7_8 + 
	0.25*om7_9)^2 + (-25.0*(ps7_10-ps7_8) * (ph8_9-ph6_9) + 25.0*(ps8_9-ps6_9) * 
	(ph7_10-ph7_8) - 200.0*ph7_9 + 100.0*ph8_9 + 100.0*ph6_9 - 200.0*ph7_9 + 
	100.0*ph7_10 + 100.0*ph7_8)^2 + (-200.0*om8_9 + 100.0*om9_9 + 100.0*om7_9 - 
	200.0*om8_9 + 100.0*om8_10 + 100.0*om8_8 + 1.5308084989341916e-13*ph9_9 - 
	1.5308084989341916e-13*ph7_9 - 2500.0*ph8_10 + 2500.0*ph8_8)^2 + (-200.0*ps8_9 
	+ 100.0*ps9_9 + 100.0*ps7_9 - 200.0*ps8_9 + 100.0*ps8_10 + 100.0*ps8_8 + 
	0.25*om8_9)^2 + (-25.0*(ps8_10-ps8_8) * (ph9_9-ph7_9) + 25.0*(ps9_9-ps7_9) * 
	(ph8_10-ph8_8) - 200.0*ph8_9 + 100.0*ph9_9 + 100.0*ph7_9 - 200.0*ph8_9 + 
	100.0*ph8_10 + 100.0*ph8_8)^2 + (-200.0*om9_9 + 100.0*om10_9 + 100.0*om8_9 - 
	200.0*om9_9 + 100.0*om9_10 + 100.0*om9_8 + 1.5308084989341916e-13*ph10_9 - 
	1.5308084989341916e-13*ph8_9 - 2500.0*ph9_10 + 2500.0*ph9_8)^2 + (-200.0*ps9_9 
	+ 100.0*ps10_9 + 100.0*ps8_9 - 200.0*ps9_9 + 100.0*ps9_10 + 100.0*ps9_8 + 
	0.25*om9_9)^2 + (-25.0*(ps9_10-ps9_8) * (ph10_9-ph8_9) + 25.0*(ps10_9-ps8_9) * 
	(ph9_10-ph9_8) - 200.0*ph9_9 + 100.0*ph10_9 + 100.0*ph8_9 - 200.0*ph9_9 + 
	100.0*ph9_10 + 100.0*ph9_8)^2 + (20.0*phm10_10 - 20.0*phm10_9 + 20.0*phm9_10 - 
	20.0*phm10_10)^2 + (20.0*phm10_9 - 20.0*phm10_10 + 20.0*phm9_10 - 
	20.0*phm10_10)^2 + (20.0*ph10_10 - 20.0*ph9_10 + 20.0*ph10_9 - 
	20.0*ph10_10 + 1.0)^2 + (-20.0*psm10_10 + 20.0*psm10_9 + 20.0*psm9_10 - 
	20.0*psm10_10)^2 + (20.0*psm10_9 - 20.0*psm10_10 + 20.0*psm9_10 - 
	20.0*psm10_10)^2 + (-20.0*ps10_10 + 20.0*ps9_10 + 20.0*ps10_9 - 
	20.0*ps10_10)^2 + (20.0*phm9_10 - 20.0*phm9_9 + 1.0)^2 + (20.0*phm9_9 - 
	20.0*phm9_10 + 1.0)^2 + (20.0*ph10_9 - 20.0*ph9_9)^2 + (20.0*phm9_9 - 
	20.0*phm10_9)^2 + (-20.0*psm9_10 + 20.0*psm9_9)^2 + (20.0*psm9_9 - 
	20.0*psm9_10)^2 + (-20.0*ps10_9 + 20.0*ps9_9)^2 + (20.0*psm9_9 - 
	20.0*psm10_9)^2 + (20.0*phm8_10 - 20.0*phm8_9 + 1.0)^2 + (20.0*phm8_9 - 
	20.0*phm8_10 + 1.0)^2 + (20.0*ph10_8 - 20.0*ph9_8)^2 + (20.0*phm9_8 - 
	20.0*phm10_8)^2 + (-20.0*psm8_10 + 20.0*psm8_9)^2 + (20.0*psm8_9 - 
	20.0*psm8_10)^2 + (-20.0*ps10_8 + 20.0*ps9_8)^2 + (20.0*psm9_8 - 
	20.0*psm10_8)^2 + (20.0*phm7_10 - 20.0*phm7_9 + 1.0)^2 + (20.0*phm7_9 - 
	20.0*phm7_10 + 1.0)^2 + (20.0*ph10_7 - 20.0*ph9_7)^2 + (20.0*phm9_7 - 
	20.0*phm10_7)^2 + (-20.0*psm7_10 + 20.0*psm7_9)^2 + (20.0*psm7_9 - 
	20.0*psm7_10)^2 + (-20.0*ps10_7 + 20.0*ps9_7)^2 + (20.0*psm9_7 - 
	20.0*psm10_7)^2 + (20.0*phm6_10 - 20.0*phm6_9 + 1.0)^2 + (20.0*phm6_9 - 
	20.0*phm6_10 + 1.0)^2 + (20.0*ph10_6 - 20.0*ph9_6)^2 + (20.0*phm9_6 - 
	20.0*phm10_6)^2 + (-20.0*psm6_10 + 20.0*psm6_9)^2 + (20.0*psm6_9 - 
	20.0*psm6_10)^2 + (-20.0*ps10_6 + 20.0*ps9_6)^2 + (20.0*psm9_6 - 
	20.0*psm10_6)^2 + (20.0*phm5_10 - 20.0*phm5_9 + 1.0)^2 + (20.0*phm5_9 - 
	20.0*phm5_10 + 1.0)^2 + (20.0*ph10_5 - 20.0*ph9_5)^2 + (20.0*phm9_5 - 
	20.0*phm10_5)^2 + (-20.0*psm5_10 + 20.0*psm5_9)^2 + (20.0*psm5_9 - 
	20.0*psm5_10)^2 + (-20.0*ps10_5 + 20.0*ps9_5)^2 + (20.0*psm9_5 - 
	20.0*psm10_5)^2 + (20.0*phm4_10 - 20.0*phm4_9 + 1.0)^2 + (20.0*phm4_9 - 
	20.0*phm4_10 + 1.0)^2 + (20.0*ph10_4 - 20.0*ph9_4)^2 + (20.0*phm9_4 - 
	20.0*phm10_4)^2 + (-20.0*psm4_10 + 20.0*psm4_9)^2 + (20.0*psm4_9 - 
	20.0*psm4_10)^2 + (-20.0*ps10_4 + 20.0*ps9_4)^2 + (20.0*psm9_4 - 
	20.0*psm10_4)^2 + (20.0*phm3_10 - 20.0*phm3_9 + 1.0)^2 + (20.0*phm3_9 - 
	20.0*phm3_10 + 1.0)^2 + (20.0*ph10_3 - 20.0*ph9_3)^2 + (20.0*phm9_3 - 
	20.0*phm10_3)^2 + (-20.0*psm3_10 + 20.0*psm3_9)^2 + (20.0*psm3_9 - 
	20.0*psm3_10)^2 + (-20.0*ps10_3 + 20.0*ps9_3)^2 + (20.0*psm9_3 - 
	20.0*psm10_3)^2 + (20.0*phm2_10 - 20.0*phm2_9 + 1.0)^2 + (20.0*phm2_9 - 
	20.0*phm2_10 + 1.0)^2 + (20.0*ph10_2 - 20.0*ph9_2)^2 + (20.0*phm9_2 - 
	20.0*phm10_2)^2 + (-20.0*psm2_10 + 20.0*psm2_9)^2 + (20.0*psm2_9 - 
	20.0*psm2_10)^2 + (-20.0*ps10_2 + 20.0*ps9_2)^2 + (20.0*psm9_2 - 
	20.0*psm10_2)^2 + (20.0*phm1_10 - 20.0*phm1_9 + 1.0)^2 + (20.0*phm1_9 - 
	20.0*phm1_10 + 1.0)^2 + (20.0*ph10_1 - 20.0*ph9_1)^2 + (20.0*phm9_1 - 
	20.0*phm10_1)^2 + (-20.0*psm1_10 + 20.0*psm1_9)^2 + (20.0*psm1_9 - 
	20.0*psm1_10)^2 + (-20.0*ps10_1 + 20.0*ps9_1)^2 + (20.0*psm9_1 - 
	20.0*psm10_1)^2 + (20.0*ph0_10 - 20.0*ph0_9 + 1.0)^2 + (20.0*ph0_9 - 
	20.0*ph0_10 + 1.0)^2 + (20.0*ph10_0 - 20.0*ph9_0)^2 + (20.0*phm9_0 - 
	20.0*phm10_0)^2 + (-20.0*ps0_10 + 20.0*ps0_9)^2 + (20.0*ps0_9 - 
	20.0*ps0_10)^2 + (-20.0*ps10_0 + 20.0*ps9_0)^2 + (20.0*psm9_0 - 
	20.0*psm10_0)^2 + (20.0*ph1_10 - 20.0*ph1_9 + 1.0)^2 + (20.0*ph1_9 - 
	20.0*ph1_10 + 1.0)^2 + (20.0*ph10_1 - 20.0*ph9_1)^2 + (20.0*phm9_1 - 
	20.0*phm10_1)^2 + (-20.0*ps1_10 + 20.0*ps1_9)^2 + (20.0*ps1_9 - 
	20.0*ps1_10)^2 + (-20.0*ps10_1 + 20.0*ps9_1)^2 + (20.0*psm9_1 - 
	20.0*psm10_1)^2 + (20.0*ph2_10 - 20.0*ph2_9 + 1.0)^2 + (20.0*ph2_9 - 
	20.0*ph2_10 + 1.0)^2 + (20.0*ph10_2 - 20.0*ph9_2)^2 + (20.0*phm9_2 - 
	20.0*phm10_2)^2 + (-20.0*ps2_10 + 20.0*ps2_9)^2 + (20.0*ps2_9 - 
	20.0*ps2_10)^2 + (-20.0*ps10_2 + 20.0*ps9_2)^2 + (20.0*psm9_2 - 
	20.0*psm10_2)^2 + (20.0*ph3_10 - 20.0*ph3_9 + 1.0)^2 + (20.0*ph3_9 - 
	20.0*ph3_10 + 1.0)^2 + (20.0*ph10_3 - 20.0*ph9_3)^2 + (20.0*phm9_3 - 
	20.0*phm10_3)^2 + (-20.0*ps3_10 + 20.0*ps3_9)^2 + (20.0*ps3_9 - 
	20.0*ps3_10)^2 + (-20.0*ps10_3 + 20.0*ps9_3)^2 + (20.0*psm9_3 - 
	20.0*psm10_3)^2 + (20.0*ph4_10 - 20.0*ph4_9 + 1.0)^2 + (20.0*ph4_9 - 
	20.0*ph4_10 + 1.0)^2 + (20.0*ph10_4 - 20.0*ph9_4)^2 + (20.0*phm9_4 - 
	20.0*phm10_4)^2 + (-20.0*ps4_10 + 20.0*ps4_9)^2 + (20.0*ps4_9 - 
	20.0*ps4_10)^2 + (-20.0*ps10_4 + 20.0*ps9_4)^2 + (20.0*psm9_4 - 
	20.0*psm10_4)^2 + (20.0*ph5_10 - 20.0*ph5_9 + 1.0)^2 + (20.0*ph5_9 - 
	20.0*ph5_10 + 1.0)^2 + (20.0*ph10_5 - 20.0*ph9_5)^2 + (20.0*phm9_5 - 
	20.0*phm10_5)^2 + (-20.0*ps5_10 + 20.0*ps5_9)^2 + (20.0*ps5_9 - 
	20.0*ps5_10)^2 + (-20.0*ps10_5 + 20.0*ps9_5)^2 + (20.0*psm9_5 - 
	20.0*psm10_5)^2 + (20.0*ph6_10 - 20.0*ph6_9 + 1.0)^2 + (20.0*ph6_9 - 
	20.0*ph6_10 + 1.0)^2 + (20.0*ph10_6 - 20.0*ph9_6)^2 + (20.0*phm9_6 - 
	20.0*phm10_6)^2 + (-20.0*ps6_10 + 20.0*ps6_9)^2 + (20.0*ps6_9 - 
	20.0*ps6_10)^2 + (-20.0*ps10_6 + 20.0*ps9_6)^2 + (20.0*psm9_6 - 
	20.0*psm10_6)^2 + (20.0*ph7_10 - 20.0*ph7_9 + 1.0)^2 + (20.0*ph7_9 - 
	20.0*ph7_10 + 1.0)^2 + (20.0*ph10_7 - 20.0*ph9_7)^2 + (20.0*phm9_7 - 
	20.0*phm10_7)^2 + (-20.0*ps7_10 + 20.0*ps7_9)^2 + (20.0*ps7_9 - 
	20.0*ps7_10)^2 + (-20.0*ps10_7 + 20.0*ps9_7)^2 + (20.0*psm9_7 - 
	20.0*psm10_7)^2 + (20.0*ph8_10 - 20.0*ph8_9 + 1.0)^2 + (20.0*ph8_9 - 
	20.0*ph8_10 + 1.0)^2 + (20.0*ph10_8 - 20.0*ph9_8)^2 + (20.0*phm9_8 - 
	20.0*phm10_8)^2 + (-20.0*ps8_10 + 20.0*ps8_9)^2 + (20.0*ps8_9 - 
	20.0*ps8_10)^2 + (-20.0*ps10_8 + 20.0*ps9_8)^2 + (20.0*psm9_8 - 
	20.0*psm10_8)^2 + (20.0*ph9_10 - 20.0*ph9_9 + 1.0)^2 + (20.0*ph9_9 - 
	20.0*ph9_10 + 1.0)^2 + (20.0*ph10_9 - 20.0*ph9_9)^2 + (20.0*phm9_9 - 
	20.0*phm10_9)^2 + (-20.0*ps9_10 + 20.0*ps9_9)^2 + (20.0*ps9_9 - 
	20.0*ps9_10)^2 + (-20.0*ps10_9 + 20.0*ps9_9)^2 + (20.0*psm9_9 - 
	20.0*psm10_9)^2 + (20.0*ph10_10 - 20.0*ph10_9 + 20.0*ph10_10 - 20.0*ph9_10)^2 + 
	(-20.0*ps10_10 + 20.0*ps10_9 - 20.0*ps10_10 + 20.0*ps9_10)^2;


solve;
	display omm10_m10;
	display phm10_m10;
	display psm10_m10;
	display omm9_m10;
	display phm9_m10;
	display psm9_m10;
	display omm8_m10;
	display phm8_m10;
	display psm8_m10;
	display omm7_m10;
	display phm7_m10;
	display psm7_m10;
	display omm6_m10;
	display phm6_m10;
	display psm6_m10;
	display omm5_m10;
	display phm5_m10;
	display psm5_m10;
	display omm4_m10;
	display phm4_m10;
	display psm4_m10;
	display omm3_m10;
	display phm3_m10;
	display psm3_m10;
	display omm2_m10;
	display phm2_m10;
	display psm2_m10;
	display omm1_m10;
	display phm1_m10;
	display psm1_m10;
	display om0_m10;
	display ph0_m10;
	display ps0_m10;
	display om1_m10;
	display ph1_m10;
	display ps1_m10;
	display om2_m10;
	display ph2_m10;
	display ps2_m10;
	display om3_m10;
	display ph3_m10;
	display ps3_m10;
	display om4_m10;
	display ph4_m10;
	display ps4_m10;
	display om5_m10;
	display ph5_m10;
	display ps5_m10;
	display om6_m10;
	display ph6_m10;
	display ps6_m10;
	display om7_m10;
	display ph7_m10;
	display ps7_m10;
	display om8_m10;
	display ph8_m10;
	display ps8_m10;
	display om9_m10;
	display ph9_m10;
	display ps9_m10;
	display om10_m10;
	display ph10_m10;
	display ps10_m10;
	display omm10_m9;
	display phm10_m9;
	display psm10_m9;
	display omm9_m9;
	display phm9_m9;
	display psm9_m9;
	display omm8_m9;
	display phm8_m9;
	display psm8_m9;
	display omm7_m9;
	display phm7_m9;
	display psm7_m9;
	display omm6_m9;
	display phm6_m9;
	display psm6_m9;
	display omm5_m9;
	display phm5_m9;
	display psm5_m9;
	display omm4_m9;
	display phm4_m9;
	display psm4_m9;
	display omm3_m9;
	display phm3_m9;
	display psm3_m9;
	display omm2_m9;
	display phm2_m9;
	display psm2_m9;
	display omm1_m9;
	display phm1_m9;
	display psm1_m9;
	display om0_m9;
	display ph0_m9;
	display ps0_m9;
	display om1_m9;
	display ph1_m9;
	display ps1_m9;
	display om2_m9;
	display ph2_m9;
	display ps2_m9;
	display om3_m9;
	display ph3_m9;
	display ps3_m9;
	display om4_m9;
	display ph4_m9;
	display ps4_m9;
	display om5_m9;
	display ph5_m9;
	display ps5_m9;
	display om6_m9;
	display ph6_m9;
	display ps6_m9;
	display om7_m9;
	display ph7_m9;
	display ps7_m9;
	display om8_m9;
	display ph8_m9;
	display ps8_m9;
	display om9_m9;
	display ph9_m9;
	display ps9_m9;
	display om10_m9;
	display ph10_m9;
	display ps10_m9;
	display omm10_m8;
	display phm10_m8;
	display psm10_m8;
	display omm9_m8;
	display phm9_m8;
	display psm9_m8;
	display omm8_m8;
	display phm8_m8;
	display psm8_m8;
	display omm7_m8;
	display phm7_m8;
	display psm7_m8;
	display omm6_m8;
	display phm6_m8;
	display psm6_m8;
	display omm5_m8;
	display phm5_m8;
	display psm5_m8;
	display omm4_m8;
	display phm4_m8;
	display psm4_m8;
	display omm3_m8;
	display phm3_m8;
	display psm3_m8;
	display omm2_m8;
	display phm2_m8;
	display psm2_m8;
	display omm1_m8;
	display phm1_m8;
	display psm1_m8;
	display om0_m8;
	display ph0_m8;
	display ps0_m8;
	display om1_m8;
	display ph1_m8;
	display ps1_m8;
	display om2_m8;
	display ph2_m8;
	display ps2_m8;
	display om3_m8;
	display ph3_m8;
	display ps3_m8;
	display om4_m8;
	display ph4_m8;
	display ps4_m8;
	display om5_m8;
	display ph5_m8;
	display ps5_m8;
	display om6_m8;
	display ph6_m8;
	display ps6_m8;
	display om7_m8;
	display ph7_m8;
	display ps7_m8;
	display om8_m8;
	display ph8_m8;
	display ps8_m8;
	display om9_m8;
	display ph9_m8;
	display ps9_m8;
	display om10_m8;
	display ph10_m8;
	display ps10_m8;
	display omm10_m7;
	display phm10_m7;
	display psm10_m7;
	display omm9_m7;
	display phm9_m7;
	display psm9_m7;
	display omm8_m7;
	display phm8_m7;
	display psm8_m7;
	display omm7_m7;
	display phm7_m7;
	display psm7_m7;
	display omm6_m7;
	display phm6_m7;
	display psm6_m7;
	display omm5_m7;
	display phm5_m7;
	display psm5_m7;
	display omm4_m7;
	display phm4_m7;
	display psm4_m7;
	display omm3_m7;
	display phm3_m7;
	display psm3_m7;
	display omm2_m7;
	display phm2_m7;
	display psm2_m7;
	display omm1_m7;
	display phm1_m7;
	display psm1_m7;
	display om0_m7;
	display ph0_m7;
	display ps0_m7;
	display om1_m7;
	display ph1_m7;
	display ps1_m7;
	display om2_m7;
	display ph2_m7;
	display ps2_m7;
	display om3_m7;
	display ph3_m7;
	display ps3_m7;
	display om4_m7;
	display ph4_m7;
	display ps4_m7;
	display om5_m7;
	display ph5_m7;
	display ps5_m7;
	display om6_m7;
	display ph6_m7;
	display ps6_m7;
	display om7_m7;
	display ph7_m7;
	display ps7_m7;
	display om8_m7;
	display ph8_m7;
	display ps8_m7;
	display om9_m7;
	display ph9_m7;
	display ps9_m7;
	display om10_m7;
	display ph10_m7;
	display ps10_m7;
	display omm10_m6;
	display phm10_m6;
	display psm10_m6;
	display omm9_m6;
	display phm9_m6;
	display psm9_m6;
	display omm8_m6;
	display phm8_m6;
	display psm8_m6;
	display omm7_m6;
	display phm7_m6;
	display psm7_m6;
	display omm6_m6;
	display phm6_m6;
	display psm6_m6;
	display omm5_m6;
	display phm5_m6;
	display psm5_m6;
	display omm4_m6;
	display phm4_m6;
	display psm4_m6;
	display omm3_m6;
	display phm3_m6;
	display psm3_m6;
	display omm2_m6;
	display phm2_m6;
	display psm2_m6;
	display omm1_m6;
	display phm1_m6;
	display psm1_m6;
	display om0_m6;
	display ph0_m6;
	display ps0_m6;
	display om1_m6;
	display ph1_m6;
	display ps1_m6;
	display om2_m6;
	display ph2_m6;
	display ps2_m6;
	display om3_m6;
	display ph3_m6;
	display ps3_m6;
	display om4_m6;
	display ph4_m6;
	display ps4_m6;
	display om5_m6;
	display ph5_m6;
	display ps5_m6;
	display om6_m6;
	display ph6_m6;
	display ps6_m6;
	display om7_m6;
	display ph7_m6;
	display ps7_m6;
	display om8_m6;
	display ph8_m6;
	display ps8_m6;
	display om9_m6;
	display ph9_m6;
	display ps9_m6;
	display om10_m6;
	display ph10_m6;
	display ps10_m6;
	display omm10_m5;
	display phm10_m5;
	display psm10_m5;
	display omm9_m5;
	display phm9_m5;
	display psm9_m5;
	display omm8_m5;
	display phm8_m5;
	display psm8_m5;
	display omm7_m5;
	display phm7_m5;
	display psm7_m5;
	display omm6_m5;
	display phm6_m5;
	display psm6_m5;
	display omm5_m5;
	display phm5_m5;
	display psm5_m5;
	display omm4_m5;
	display phm4_m5;
	display psm4_m5;
	display omm3_m5;
	display phm3_m5;
	display psm3_m5;
	display omm2_m5;
	display phm2_m5;
	display psm2_m5;
	display omm1_m5;
	display phm1_m5;
	display psm1_m5;
	display om0_m5;
	display ph0_m5;
	display ps0_m5;
	display om1_m5;
	display ph1_m5;
	display ps1_m5;
	display om2_m5;
	display ph2_m5;
	display ps2_m5;
	display om3_m5;
	display ph3_m5;
	display ps3_m5;
	display om4_m5;
	display ph4_m5;
	display ps4_m5;
	display om5_m5;
	display ph5_m5;
	display ps5_m5;
	display om6_m5;
	display ph6_m5;
	display ps6_m5;
	display om7_m5;
	display ph7_m5;
	display ps7_m5;
	display om8_m5;
	display ph8_m5;
	display ps8_m5;
	display om9_m5;
	display ph9_m5;
	display ps9_m5;
	display om10_m5;
	display ph10_m5;
	display ps10_m5;
	display omm10_m4;
	display phm10_m4;
	display psm10_m4;
	display omm9_m4;
	display phm9_m4;
	display psm9_m4;
	display omm8_m4;
	display phm8_m4;
	display psm8_m4;
	display omm7_m4;
	display phm7_m4;
	display psm7_m4;
	display omm6_m4;
	display phm6_m4;
	display psm6_m4;
	display omm5_m4;
	display phm5_m4;
	display psm5_m4;
	display omm4_m4;
	display phm4_m4;
	display psm4_m4;
	display omm3_m4;
	display phm3_m4;
	display psm3_m4;
	display omm2_m4;
	display phm2_m4;
	display psm2_m4;
	display omm1_m4;
	display phm1_m4;
	display psm1_m4;
	display om0_m4;
	display ph0_m4;
	display ps0_m4;
	display om1_m4;
	display ph1_m4;
	display ps1_m4;
	display om2_m4;
	display ph2_m4;
	display ps2_m4;
	display om3_m4;
	display ph3_m4;
	display ps3_m4;
	display om4_m4;
	display ph4_m4;
	display ps4_m4;
	display om5_m4;
	display ph5_m4;
	display ps5_m4;
	display om6_m4;
	display ph6_m4;
	display ps6_m4;
	display om7_m4;
	display ph7_m4;
	display ps7_m4;
	display om8_m4;
	display ph8_m4;
	display ps8_m4;
	display om9_m4;
	display ph9_m4;
	display ps9_m4;
	display om10_m4;
	display ph10_m4;
	display ps10_m4;
	display omm10_m3;
	display phm10_m3;
	display psm10_m3;
	display omm9_m3;
	display phm9_m3;
	display psm9_m3;
	display omm8_m3;
	display phm8_m3;
	display psm8_m3;
	display omm7_m3;
	display phm7_m3;
	display psm7_m3;
	display omm6_m3;
	display phm6_m3;
	display psm6_m3;
	display omm5_m3;
	display phm5_m3;
	display psm5_m3;
	display omm4_m3;
	display phm4_m3;
	display psm4_m3;
	display omm3_m3;
	display phm3_m3;
	display psm3_m3;
	display omm2_m3;
	display phm2_m3;
	display psm2_m3;
	display omm1_m3;
	display phm1_m3;
	display psm1_m3;
	display om0_m3;
	display ph0_m3;
	display ps0_m3;
	display om1_m3;
	display ph1_m3;
	display ps1_m3;
	display om2_m3;
	display ph2_m3;
	display ps2_m3;
	display om3_m3;
	display ph3_m3;
	display ps3_m3;
	display om4_m3;
	display ph4_m3;
	display ps4_m3;
	display om5_m3;
	display ph5_m3;
	display ps5_m3;
	display om6_m3;
	display ph6_m3;
	display ps6_m3;
	display om7_m3;
	display ph7_m3;
	display ps7_m3;
	display om8_m3;
	display ph8_m3;
	display ps8_m3;
	display om9_m3;
	display ph9_m3;
	display ps9_m3;
	display om10_m3;
	display ph10_m3;
	display ps10_m3;
	display omm10_m2;
	display phm10_m2;
	display psm10_m2;
	display omm9_m2;
	display phm9_m2;
	display psm9_m2;
	display omm8_m2;
	display phm8_m2;
	display psm8_m2;
	display omm7_m2;
	display phm7_m2;
	display psm7_m2;
	display omm6_m2;
	display phm6_m2;
	display psm6_m2;
	display omm5_m2;
	display phm5_m2;
	display psm5_m2;
	display omm4_m2;
	display phm4_m2;
	display psm4_m2;
	display omm3_m2;
	display phm3_m2;
	display psm3_m2;
	display omm2_m2;
	display phm2_m2;
	display psm2_m2;
	display omm1_m2;
	display phm1_m2;
	display psm1_m2;
	display om0_m2;
	display ph0_m2;
	display ps0_m2;
	display om1_m2;
	display ph1_m2;
	display ps1_m2;
	display om2_m2;
	display ph2_m2;
	display ps2_m2;
	display om3_m2;
	display ph3_m2;
	display ps3_m2;
	display om4_m2;
	display ph4_m2;
	display ps4_m2;
	display om5_m2;
	display ph5_m2;
	display ps5_m2;
	display om6_m2;
	display ph6_m2;
	display ps6_m2;
	display om7_m2;
	display ph7_m2;
	display ps7_m2;
	display om8_m2;
	display ph8_m2;
	display ps8_m2;
	display om9_m2;
	display ph9_m2;
	display ps9_m2;
	display om10_m2;
	display ph10_m2;
	display ps10_m2;
	display omm10_m1;
	display phm10_m1;
	display psm10_m1;
	display omm9_m1;
	display phm9_m1;
	display psm9_m1;
	display omm8_m1;
	display phm8_m1;
	display psm8_m1;
	display omm7_m1;
	display phm7_m1;
	display psm7_m1;
	display omm6_m1;
	display phm6_m1;
	display psm6_m1;
	display omm5_m1;
	display phm5_m1;
	display psm5_m1;
	display omm4_m1;
	display phm4_m1;
	display psm4_m1;
	display omm3_m1;
	display phm3_m1;
	display psm3_m1;
	display omm2_m1;
	display phm2_m1;
	display psm2_m1;
	display omm1_m1;
	display phm1_m1;
	display psm1_m1;
	display om0_m1;
	display ph0_m1;
	display ps0_m1;
	display om1_m1;
	display ph1_m1;
	display ps1_m1;
	display om2_m1;
	display ph2_m1;
	display ps2_m1;
	display om3_m1;
	display ph3_m1;
	display ps3_m1;
	display om4_m1;
	display ph4_m1;
	display ps4_m1;
	display om5_m1;
	display ph5_m1;
	display ps5_m1;
	display om6_m1;
	display ph6_m1;
	display ps6_m1;
	display om7_m1;
	display ph7_m1;
	display ps7_m1;
	display om8_m1;
	display ph8_m1;
	display ps8_m1;
	display om9_m1;
	display ph9_m1;
	display ps9_m1;
	display om10_m1;
	display ph10_m1;
	display ps10_m1;
	display omm10_0;
	display phm10_0;
	display psm10_0;
	display omm9_0;
	display phm9_0;
	display psm9_0;
	display omm8_0;
	display phm8_0;
	display psm8_0;
	display omm7_0;
	display phm7_0;
	display psm7_0;
	display omm6_0;
	display phm6_0;
	display psm6_0;
	display omm5_0;
	display phm5_0;
	display psm5_0;
	display omm4_0;
	display phm4_0;
	display psm4_0;
	display omm3_0;
	display phm3_0;
	display psm3_0;
	display omm2_0;
	display phm2_0;
	display psm2_0;
	display omm1_0;
	display phm1_0;
	display psm1_0;
	display om0_0;
	display ph0_0;
	display ps0_0;
	display om1_0;
	display ph1_0;
	display ps1_0;
	display om2_0;
	display ph2_0;
	display ps2_0;
	display om3_0;
	display ph3_0;
	display ps3_0;
	display om4_0;
	display ph4_0;
	display ps4_0;
	display om5_0;
	display ph5_0;
	display ps5_0;
	display om6_0;
	display ph6_0;
	display ps6_0;
	display om7_0;
	display ph7_0;
	display ps7_0;
	display om8_0;
	display ph8_0;
	display ps8_0;
	display om9_0;
	display ph9_0;
	display ps9_0;
	display om10_0;
	display ph10_0;
	display ps10_0;
	display omm10_1;
	display phm10_1;
	display psm10_1;
	display omm9_1;
	display phm9_1;
	display psm9_1;
	display omm8_1;
	display phm8_1;
	display psm8_1;
	display omm7_1;
	display phm7_1;
	display psm7_1;
	display omm6_1;
	display phm6_1;
	display psm6_1;
	display omm5_1;
	display phm5_1;
	display psm5_1;
	display omm4_1;
	display phm4_1;
	display psm4_1;
	display omm3_1;
	display phm3_1;
	display psm3_1;
	display omm2_1;
	display phm2_1;
	display psm2_1;
	display omm1_1;
	display phm1_1;
	display psm1_1;
	display om0_1;
	display ph0_1;
	display ps0_1;
	display om1_1;
	display ph1_1;
	display ps1_1;
	display om2_1;
	display ph2_1;
	display ps2_1;
	display om3_1;
	display ph3_1;
	display ps3_1;
	display om4_1;
	display ph4_1;
	display ps4_1;
	display om5_1;
	display ph5_1;
	display ps5_1;
	display om6_1;
	display ph6_1;
	display ps6_1;
	display om7_1;
	display ph7_1;
	display ps7_1;
	display om8_1;
	display ph8_1;
	display ps8_1;
	display om9_1;
	display ph9_1;
	display ps9_1;
	display om10_1;
	display ph10_1;
	display ps10_1;
	display omm10_2;
	display phm10_2;
	display psm10_2;
	display omm9_2;
	display phm9_2;
	display psm9_2;
	display omm8_2;
	display phm8_2;
	display psm8_2;
	display omm7_2;
	display phm7_2;
	display psm7_2;
	display omm6_2;
	display phm6_2;
	display psm6_2;
	display omm5_2;
	display phm5_2;
	display psm5_2;
	display omm4_2;
	display phm4_2;
	display psm4_2;
	display omm3_2;
	display phm3_2;
	display psm3_2;
	display omm2_2;
	display phm2_2;
	display psm2_2;
	display omm1_2;
	display phm1_2;
	display psm1_2;
	display om0_2;
	display ph0_2;
	display ps0_2;
	display om1_2;
	display ph1_2;
	display ps1_2;
	display om2_2;
	display ph2_2;
	display ps2_2;
	display om3_2;
	display ph3_2;
	display ps3_2;
	display om4_2;
	display ph4_2;
	display ps4_2;
	display om5_2;
	display ph5_2;
	display ps5_2;
	display om6_2;
	display ph6_2;
	display ps6_2;
	display om7_2;
	display ph7_2;
	display ps7_2;
	display om8_2;
	display ph8_2;
	display ps8_2;
	display om9_2;
	display ph9_2;
	display ps9_2;
	display om10_2;
	display ph10_2;
	display ps10_2;
	display omm10_3;
	display phm10_3;
	display psm10_3;
	display omm9_3;
	display phm9_3;
	display psm9_3;
	display omm8_3;
	display phm8_3;
	display psm8_3;
	display omm7_3;
	display phm7_3;
	display psm7_3;
	display omm6_3;
	display phm6_3;
	display psm6_3;
	display omm5_3;
	display phm5_3;
	display psm5_3;
	display omm4_3;
	display phm4_3;
	display psm4_3;
	display omm3_3;
	display phm3_3;
	display psm3_3;
	display omm2_3;
	display phm2_3;
	display psm2_3;
	display omm1_3;
	display phm1_3;
	display psm1_3;
	display om0_3;
	display ph0_3;
	display ps0_3;
	display om1_3;
	display ph1_3;
	display ps1_3;
	display om2_3;
	display ph2_3;
	display ps2_3;
	display om3_3;
	display ph3_3;
	display ps3_3;
	display om4_3;
	display ph4_3;
	display ps4_3;
	display om5_3;
	display ph5_3;
	display ps5_3;
	display om6_3;
	display ph6_3;
	display ps6_3;
	display om7_3;
	display ph7_3;
	display ps7_3;
	display om8_3;
	display ph8_3;
	display ps8_3;
	display om9_3;
	display ph9_3;
	display ps9_3;
	display om10_3;
	display ph10_3;
	display ps10_3;
	display omm10_4;
	display phm10_4;
	display psm10_4;
	display omm9_4;
	display phm9_4;
	display psm9_4;
	display omm8_4;
	display phm8_4;
	display psm8_4;
	display omm7_4;
	display phm7_4;
	display psm7_4;
	display omm6_4;
	display phm6_4;
	display psm6_4;
	display omm5_4;
	display phm5_4;
	display psm5_4;
	display omm4_4;
	display phm4_4;
	display psm4_4;
	display omm3_4;
	display phm3_4;
	display psm3_4;
	display omm2_4;
	display phm2_4;
	display psm2_4;
	display omm1_4;
	display phm1_4;
	display psm1_4;
	display om0_4;
	display ph0_4;
	display ps0_4;
	display om1_4;
	display ph1_4;
	display ps1_4;
	display om2_4;
	display ph2_4;
	display ps2_4;
	display om3_4;
	display ph3_4;
	display ps3_4;
	display om4_4;
	display ph4_4;
	display ps4_4;
	display om5_4;
	display ph5_4;
	display ps5_4;
	display om6_4;
	display ph6_4;
	display ps6_4;
	display om7_4;
	display ph7_4;
	display ps7_4;
	display om8_4;
	display ph8_4;
	display ps8_4;
	display om9_4;
	display ph9_4;
	display ps9_4;
	display om10_4;
	display ph10_4;
	display ps10_4;
	display omm10_5;
	display phm10_5;
	display psm10_5;
	display omm9_5;
	display phm9_5;
	display psm9_5;
	display omm8_5;
	display phm8_5;
	display psm8_5;
	display omm7_5;
	display phm7_5;
	display psm7_5;
	display omm6_5;
	display phm6_5;
	display psm6_5;
	display omm5_5;
	display phm5_5;
	display psm5_5;
	display omm4_5;
	display phm4_5;
	display psm4_5;
	display omm3_5;
	display phm3_5;
	display psm3_5;
	display omm2_5;
	display phm2_5;
	display psm2_5;
	display omm1_5;
	display phm1_5;
	display psm1_5;
	display om0_5;
	display ph0_5;
	display ps0_5;
	display om1_5;
	display ph1_5;
	display ps1_5;
	display om2_5;
	display ph2_5;
	display ps2_5;
	display om3_5;
	display ph3_5;
	display ps3_5;
	display om4_5;
	display ph4_5;
	display ps4_5;
	display om5_5;
	display ph5_5;
	display ps5_5;
	display om6_5;
	display ph6_5;
	display ps6_5;
	display om7_5;
	display ph7_5;
	display ps7_5;
	display om8_5;
	display ph8_5;
	display ps8_5;
	display om9_5;
	display ph9_5;
	display ps9_5;
	display om10_5;
	display ph10_5;
	display ps10_5;
	display omm10_6;
	display phm10_6;
	display psm10_6;
	display omm9_6;
	display phm9_6;
	display psm9_6;
	display omm8_6;
	display phm8_6;
	display psm8_6;
	display omm7_6;
	display phm7_6;
	display psm7_6;
	display omm6_6;
	display phm6_6;
	display psm6_6;
	display omm5_6;
	display phm5_6;
	display psm5_6;
	display omm4_6;
	display phm4_6;
	display psm4_6;
	display omm3_6;
	display phm3_6;
	display psm3_6;
	display omm2_6;
	display phm2_6;
	display psm2_6;
	display omm1_6;
	display phm1_6;
	display psm1_6;
	display om0_6;
	display ph0_6;
	display ps0_6;
	display om1_6;
	display ph1_6;
	display ps1_6;
	display om2_6;
	display ph2_6;
	display ps2_6;
	display om3_6;
	display ph3_6;
	display ps3_6;
	display om4_6;
	display ph4_6;
	display ps4_6;
	display om5_6;
	display ph5_6;
	display ps5_6;
	display om6_6;
	display ph6_6;
	display ps6_6;
	display om7_6;
	display ph7_6;
	display ps7_6;
	display om8_6;
	display ph8_6;
	display ps8_6;
	display om9_6;
	display ph9_6;
	display ps9_6;
	display om10_6;
	display ph10_6;
	display ps10_6;
	display omm10_7;
	display phm10_7;
	display psm10_7;
	display omm9_7;
	display phm9_7;
	display psm9_7;
	display omm8_7;
	display phm8_7;
	display psm8_7;
	display omm7_7;
	display phm7_7;
	display psm7_7;
	display omm6_7;
	display phm6_7;
	display psm6_7;
	display omm5_7;
	display phm5_7;
	display psm5_7;
	display omm4_7;
	display phm4_7;
	display psm4_7;
	display omm3_7;
	display phm3_7;
	display psm3_7;
	display omm2_7;
	display phm2_7;
	display psm2_7;
	display omm1_7;
	display phm1_7;
	display psm1_7;
	display om0_7;
	display ph0_7;
	display ps0_7;
	display om1_7;
	display ph1_7;
	display ps1_7;
	display om2_7;
	display ph2_7;
	display ps2_7;
	display om3_7;
	display ph3_7;
	display ps3_7;
	display om4_7;
	display ph4_7;
	display ps4_7;
	display om5_7;
	display ph5_7;
	display ps5_7;
	display om6_7;
	display ph6_7;
	display ps6_7;
	display om7_7;
	display ph7_7;
	display ps7_7;
	display om8_7;
	display ph8_7;
	display ps8_7;
	display om9_7;
	display ph9_7;
	display ps9_7;
	display om10_7;
	display ph10_7;
	display ps10_7;
	display omm10_8;
	display phm10_8;
	display psm10_8;
	display omm9_8;
	display phm9_8;
	display psm9_8;
	display omm8_8;
	display phm8_8;
	display psm8_8;
	display omm7_8;
	display phm7_8;
	display psm7_8;
	display omm6_8;
	display phm6_8;
	display psm6_8;
	display omm5_8;
	display phm5_8;
	display psm5_8;
	display omm4_8;
	display phm4_8;
	display psm4_8;
	display omm3_8;
	display phm3_8;
	display psm3_8;
	display omm2_8;
	display phm2_8;
	display psm2_8;
	display omm1_8;
	display phm1_8;
	display psm1_8;
	display om0_8;
	display ph0_8;
	display ps0_8;
	display om1_8;
	display ph1_8;
	display ps1_8;
	display om2_8;
	display ph2_8;
	display ps2_8;
	display om3_8;
	display ph3_8;
	display ps3_8;
	display om4_8;
	display ph4_8;
	display ps4_8;
	display om5_8;
	display ph5_8;
	display ps5_8;
	display om6_8;
	display ph6_8;
	display ps6_8;
	display om7_8;
	display ph7_8;
	display ps7_8;
	display om8_8;
	display ph8_8;
	display ps8_8;
	display om9_8;
	display ph9_8;
	display ps9_8;
	display om10_8;
	display ph10_8;
	display ps10_8;
	display omm10_9;
	display phm10_9;
	display psm10_9;
	display omm9_9;
	display phm9_9;
	display psm9_9;
	display omm8_9;
	display phm8_9;
	display psm8_9;
	display omm7_9;
	display phm7_9;
	display psm7_9;
	display omm6_9;
	display phm6_9;
	display psm6_9;
	display omm5_9;
	display phm5_9;
	display psm5_9;
	display omm4_9;
	display phm4_9;
	display psm4_9;
	display omm3_9;
	display phm3_9;
	display psm3_9;
	display omm2_9;
	display phm2_9;
	display psm2_9;
	display omm1_9;
	display phm1_9;
	display psm1_9;
	display om0_9;
	display ph0_9;
	display ps0_9;
	display om1_9;
	display ph1_9;
	display ps1_9;
	display om2_9;
	display ph2_9;
	display ps2_9;
	display om3_9;
	display ph3_9;
	display ps3_9;
	display om4_9;
	display ph4_9;
	display ps4_9;
	display om5_9;
	display ph5_9;
	display ps5_9;
	display om6_9;
	display ph6_9;
	display ps6_9;
	display om7_9;
	display ph7_9;
	display ps7_9;
	display om8_9;
	display ph8_9;
	display ps8_9;
	display om9_9;
	display ph9_9;
	display ps9_9;
	display om10_9;
	display ph10_9;
	display ps10_9;
	display omm10_10;
	display phm10_10;
	display psm10_10;
	display omm9_10;
	display phm9_10;
	display psm9_10;
	display omm8_10;
	display phm8_10;
	display psm8_10;
	display omm7_10;
	display phm7_10;
	display psm7_10;
	display omm6_10;
	display phm6_10;
	display psm6_10;
	display omm5_10;
	display phm5_10;
	display psm5_10;
	display omm4_10;
	display phm4_10;
	display psm4_10;
	display omm3_10;
	display phm3_10;
	display psm3_10;
	display omm2_10;
	display phm2_10;
	display psm2_10;
	display omm1_10;
	display phm1_10;
	display psm1_10;
	display om0_10;
	display ph0_10;
	display ps0_10;
	display om1_10;
	display ph1_10;
	display ps1_10;
	display om2_10;
	display ph2_10;
	display ps2_10;
	display om3_10;
	display ph3_10;
	display ps3_10;
	display om4_10;
	display ph4_10;
	display ps4_10;
	display om5_10;
	display ph5_10;
	display ps5_10;
	display om6_10;
	display ph6_10;
	display ps6_10;
	display om7_10;
	display ph7_10;
	display ps7_10;
	display om8_10;
	display ph8_10;
	display ps8_10;
	display om9_10;
	display ph9_10;
	display ps9_10;
	display om10_10;
	display ph10_10;
	display ps10_10;
display obj;
