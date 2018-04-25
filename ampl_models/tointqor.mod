#
#**************************
# SET UP THE INITIAL DATA *
#**************************
#   Problem :
#   *********
#   Toint's  Quadratic Operations Research problem
#   Source:
#   Ph.L. Toint,
#   "Some numerical results using a sparse matrix updating formula in
#   unconstrained optimization",
#   Mathematics of Computation 32(1):839-852, 1978.
#   See also Buckley#55 (p.94) (With a slightly lower optimal value?)
#   SIF input: Ph. Toint, Dec 1989.
#   classification QUR2-MN-50-0
#   Number of variables
#   Problem parameters
#   Problem data
#   Other parameters
#   Solution
	param n := 50;
	param alph1 := 1.25;
	param alph2 := 1.4;
	param alph3 := 2.4;
	param alph4 := 1.4;
	param alph5 := 1.75;
	param alph6 := 1.2;
	param alph7 := 2.25;
	param alph8 := 1.2;
	param alph9 := 1.0;
	param alph10 := 1.1;
	param alph11 := 1.5;
	param alph12 := 1.6;
	param alph13 := 1.25;
	param alph14 := 1.25;
	param alph15 := 1.2;
	param alph16 := 1.2;
	param alph17 := 1.4;
	param alph18 := 0.5;
	param alph19 := 0.5;
	param alph20 := 1.25;
	param alph21 := 1.8;
	param alph22 := 0.75;
	param alph23 := 1.25;
	param alph24 := 1.4;
	param alph25 := 1.6;
	param alph26 := 2.0;
	param alph27 := 1.0;
	param alph28 := 1.6;
	param alph29 := 1.25;
	param alph30 := 2.75;
	param alph31 := 1.25;
	param alph32 := 1.25;
	param alph33 := 1.25;
	param alph34 := 3.0;
	param alph35 := 1.5;
	param alph36 := 2.0;
	param alph37 := 1.25;
	param alph38 := 1.4;
	param alph39 := 1.8;
	param alph40 := 1.5;
	param alph41 := 2.2;
	param alph42 := 1.4;
	param alph43 := 1.5;
	param alph44 := 1.25;
	param alph45 := 2.0;
	param alph46 := 1.5;
	param alph47 := 1.25;
	param alph48 := 1.4;
	param alph49 := 0.6;
	param alph50 := 1.5;
	param beta1 := 1.0;
	param beta2 := 1.5;
	param beta3 := 1.0;
	param beta4 := 0.1;
	param beta5 := 1.5;
	param beta6 := 2.0;
	param beta7 := 1.0;
	param beta8 := 1.5;
	param beta9 := 3.0;
	param beta10 := 2.0;
	param beta11 := 1.0;
	param beta12 := 3.0;
	param beta13 := 0.1;
	param beta14 := 1.5;
	param beta15 := 0.15;
	param beta16 := 2.0;
	param beta17 := 1.0;
	param beta18 := 0.1;
	param beta19 := 3.0;
	param beta20 := 0.1;
	param beta21 := 1.2;
	param beta22 := 1.0;
	param beta23 := 0.1;
	param beta24 := 2.0;
	param beta25 := 1.2;
	param beta26 := 3.0;
	param beta27 := 1.5;
	param beta28 := 3.0;
	param beta29 := 2.0;
	param beta30 := 1.0;
	param beta31 := 1.2;
	param beta32 := 2.0;
	param beta33 := 1.0;
	param d1 := -5.0;
	param d2 := -5.0;
	param d3 := -5.0;
	param d4 := -2.5;
	param d5 := -6.0;
	param d6 := -6.0;
	param d7 := -5.0;
	param d8 := -6.0;
	param d9 := -10.0;
	param d10 := -6.0;
	param d11 := -5.0;
	param d12 := -9.0;
	param d13 := -2.0;
	param d14 := -7.0;
	param d15 := -2.5;
	param d16 := -6.0;
	param d17 := -5.0;
	param d18 := -2.0;
	param d19 := -9.0;
	param d20 := -2.0;
	param d21 := -5.0;
	param d22 := -5.0;
	param d23 := -2.5;
	param d24 := -5.0;
	param d25 := -6.0;
	param d26 := -10.0;
	param d27 := -7.0;
	param d28 := -10.0;
	param d29 := -6.0;
	param d30 := -5.0;
	param d31 := -4.0;
	param d32 := -4.0;
	param d33 := -4.0;
	param scale := 1.0 / (1.0);

	var x1;
	var x2;
	var x3;
	var x4;
	var x5;
	var x6;
	var x7;
	var x8;
	var x9;
	var x10;
	var x11;
	var x12;
	var x13;
	var x14;
	var x15;
	var x16;
	var x17;
	var x18;
	var x19;
	var x20;
	var x21;
	var x22;
	var x23;
	var x24;
	var x25;
	var x26;
	var x27;
	var x28;
	var x29;
	var x30;
	var x31;
	var x32;
	var x33;
	var x34;
	var x35;
	var x36;
	var x37;
	var x38;
	var x39;
	var x40;
	var x41;
	var x42;
	var x43;
	var x44;
	var x45;
	var x46;
	var x47;
	var x48;
	var x49;
	var x50;

minimize obj:
	(((x1)*(x1))/0.8) + (((x2)*(x2))/0.7142857142857143) + 
	(((x3)*(x3))/0.4166666666666667) + (((x4)*(x4))/0.7142857142857143) + 
	(((x5)*(x5))/0.5714285714285714) + (((x6)*(x6))/0.8333333333333334) + 
	(((x7)*(x7))/0.4444444444444444) + (((x8)*(x8))/0.8333333333333334) + (x9)*(x9) 
	+ (((x10)*(x10))/0.9090909090909091) + (((x11)*(x11))/0.6666666666666666) + 
	(((x12)*(x12))/0.625) + (((x13)*(x13))/0.8) + (((x14)*(x14))/0.8) + 
	(((x15)*(x15))/0.8333333333333334) + (((x16)*(x16))/0.8333333333333334) + 
	(((x17)*(x17))/0.7142857142857143) + (((x18)*(x18))/2.0) + (((x19)*(x19))/2.0) 
	+ (((x20)*(x20))/0.8) + (((x21)*(x21))/0.5555555555555556) + 
	(((x22)*(x22))/1.3333333333333333) + (((x23)*(x23))/0.8) + 
	(((x24)*(x24))/0.7142857142857143) + (((x25)*(x25))/0.625) + 
	(((x26)*(x26))/0.5) + (x27)*(x27) + (((x28)*(x28))/0.625) + (((x29)*(x29))/0.8) 
	+ (((x30)*(x30))/0.36363636363636365) + (((x31)*(x31))/0.8) + 
	(((x32)*(x32))/0.8) + (((x33)*(x33))/0.8) + (((x34)*(x34))/0.3333333333333333) 
	+ (((x35)*(x35))/0.6666666666666666) + (((x36)*(x36))/0.5) + 
	(((x37)*(x37))/0.8) + (((x38)*(x38))/0.7142857142857143) + 
	(((x39)*(x39))/0.5555555555555556) + (((x40)*(x40))/0.6666666666666666) + 
	(((x41)*(x41))/0.45454545454545453) + (((x42)*(x42))/0.7142857142857143) + 
	(((x43)*(x43))/0.6666666666666666) + (((x44)*(x44))/0.8) + (((x45)*(x45))/0.5) 
	+ (((x46)*(x46))/0.6666666666666666) + (((x47)*(x47))/0.8) + 
	(((x48)*(x48))/0.7142857142857143) + (((x49)*(x49))/1.6666666666666667) + 
	(((x50)*(x50))/0.6666666666666666) + (-x31 + x1 + 5.0)*(-x31 + x1 + 5.0) + 
	(((-x1 + x2 + x3 + 5.0)*(-x1 + x2 + x3 + 5.0))/0.6666666666666666) + (-x2 + x4 
	+ x5 + 5.0)*(-x2 + x4 + x5 + 5.0) + (((-x4 + x6 + x7 + 2.5)*(-x4 + x6 + x7 + 
	2.5))/10.0) + (((-x6 + x8 + x9 + 6.0)*(-x6 + x8 + x9 + 
	6.0))/0.6666666666666666) + (((-x8 + x10 + x11 + 6.0)*(-x8 + x10 + x11 + 
	6.0))/0.5) + (-x10 + x12 + x13 + 5.0)*(-x10 + x12 + x13 + 5.0) + (((-x12 + x14 
	+ x15 + 6.0)*(-x12 + x14 + x15 + 6.0))/0.6666666666666666) + (((-x11 - x13 - 
	x14 + x16 + x17 + 10.0)*(-x11 - x13 - x14 + x16 + x17 + 
	10.0))/0.3333333333333333) + (((-x16 + x18 + x19 + 6.0)*(-x16 + x18 + x19 + 
	6.0))/0.5) + (-x9 - x18 + x20 + 5.0)*(-x9 - x18 + x20 + 5.0) + (((-x5 - x20 - 
	x21 + 9.0)*(-x5 - x20 - x21 + 9.0))/0.3333333333333333) + (((-x19 + x22 + x23 + 
	x24 + 2.0)*(-x19 + x22 + x23 + x24 + 2.0))/10.0) + (((-x23 + x25 + x26 + 
	7.0)*(-x23 + x25 + x26 + 7.0))/0.6666666666666666) + (((-x7 - x25 + x27 + x28 + 
	2.5)*(-x7 - x25 + x27 + x28 + 2.5))/6.666666666666667) + (((-x28 + x29 + x30 + 
	6.0)*(-x28 + x29 + x30 + 6.0))/0.5) + (-x29 + x31 + x32 + 5.0)*(-x29 + x31 + 
	x32 + 5.0) + (((-x32 + x33 + x34 + 2.0)*(-x32 + x33 + x34 + 2.0))/10.0) + 
	(((-x3 - x33 + x35 + 9.0)*(-x3 - x33 + x35 + 9.0))/0.3333333333333333) + 
	(((-x35 + x21 + x36 + 2.0)*(-x35 + x21 + x36 + 2.0))/10.0) + (((-x36 + x37 + 
	x38 + 5.0)*(-x36 + x37 + x38 + 5.0))/0.8333333333333334) + (-x30 - x37 + x39 + 
	5.0)*(-x30 - x37 + x39 + 5.0) + (((-x38 - x39 + x40 + 2.5)*(-x38 - x39 + x40 + 
	2.5))/10.0) + (((-x40 + x41 + x42 + 5.0)*(-x40 + x41 + x42 + 5.0))/0.5) + 
	(((-x41 + x43 + x44 + x50 + 6.0)*(-x41 + x43 + x44 + x50 + 
	6.0))/0.8333333333333334) + (((-x44 + x45 + x46 + x47 + 10.0)*(-x44 + x45 + x46 
	+ x47 + 10.0))/0.3333333333333333) + (((-x46 + x48 + 7.0)*(-x46 + x48 + 
	7.0))/0.6666666666666666) + (((-x42 - x45 - x48 - x50 + x49 + 10.0)*(-x42 - x45 
	- x48 - x50 + x49 + 10.0))/0.3333333333333333) + (((-x26 - x34 - x43 + 
	6.0)*(-x26 - x34 - x43 + 6.0))/0.5) + (-x15 - x17 - x24 - x47 + 5.0)*(-x15 - 
	x17 - x24 - x47 + 5.0) + (((-x49 + 4.0)*(-x49 + 4.0))/0.8333333333333334) + 
	(((-x22 + 4.0)*(-x22 + 4.0))/0.5) + (-x27 + 4.0)*(-x27 + 4.0);


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
	display x10;
	display x11;
	display x12;
	display x13;
	display x14;
	display x15;
	display x16;
	display x17;
	display x18;
	display x19;
	display x20;
	display x21;
	display x22;
	display x23;
	display x24;
	display x25;
	display x26;
	display x27;
	display x28;
	display x29;
	display x30;
	display x31;
	display x32;
	display x33;
	display x34;
	display x35;
	display x36;
	display x37;
	display x38;
	display x39;
	display x40;
	display x41;
	display x42;
	display x43;
	display x44;
	display x45;
	display x46;
	display x47;
	display x48;
	display x49;
	display x50;
display obj;
