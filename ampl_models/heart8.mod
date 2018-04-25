#
#**************************
# SET UP THE INITIAL DATA *
#**************************
#   Problem :
#   *********
#   Dipole model of the heart (6 x 6 version).
#   Source:
#   J. E. Dennis, Jr., D. M. Gay, P. A. Vu,
#   "A New Nonlinear Equations Test Problem".
#   Tech. Rep. 83-16, Dept. of Math. Sci., Rice Univ., Houston, TX
#   June 1983, revised May 1985.
#   SIF input: A.R. Conn, May 1993.
#              correction by Ph. Shott, January, 1995.
#   classification NOR2-MN-8-8
#   some useful parameters.
#RE sum_Mx              0.485
#RE sum_My              -0.0019
#RE sum_A               -0.0581
#RE sum_B               0.015
#RE sum_C               0.105
#RE sum_D               0.0406
#RE sum_E               0.167
#RE sum_F               -0.399
#RE sum_Mx              -0.816
#RE sum_My              -0.017
#RE sum_A               -1.826
#RE sum_B               -0.754
#RE sum_C               -4.839
#RE sum_D               -3.259
#RE sum_E               -14.023
#RE sum_F               15.467
#RE sum_Mx              -0.809
#RE sum_My              -0.021
#RE sum_A               -2.04
#RE sum_B               -0.614
#RE sum_C               -6.903
#RE sum_D               -2.934
#RE sum_E               -26.328
#RE sum_F               18.639
#RE sum_Mx              -0.807
#RE sum_My              -0.021
#RE sum_A               -2.379
#RE sum_B               -0.364
#RE sum_C               -10.541
#RE sum_D               -1.961
#RE sum_E               -51.551
#RE sum_F               21.053
#   Solution
	param sum_mx := -0.69;
	param sum_my := -0.044;
	param sum_a := -1.57;
	param sum_b := -1.31;
	param sum_c := -2.65;
	param sum_d := 2.0;
	param sum_e := -12.6;
	param sum_f := 9.48;

	var a;
	var b := 1.0;
	var c;
	var d := 1.0;
	var t := 1.0;
	var u := 1.0;
	var v := 1.0;
	var w := 1.0;

minimize obj: 0;
subject to cons1:
	(a + b + 0.69) = 0;
subject to cons2:
	(c + d + 0.044) = 0;
subject to cons3:
	(t * a + u * b - v * c - w * d + 1.57) = 0;
subject to cons4:
	(v * a + w * b + t * c + u * d + 1.31) = 0;
subject to cons5:
	(a * (t^2-v^2) - 2.0*c * t * v + b * (u^2-w^2) - 2.0*d * u * w + 2.65) = 0;
subject to cons6:
	(c * (t^2-v^2) + 2.0*a * t * v + d * (u^2-w^2) + 2.0*b * u * w - 2.0) = 0;
subject to cons7:
	(a * t * (t^2-(3.0)*v^2) + c * v * (v^2-(3.0)*t^2) + b * u * (u^2-(3.0)*w^2) + d * w * (w^2-(3.0)*u^2) + 12.6) = 0;
subject to cons8:
	(c * t * (t^2-(3.0)*v^2) - a * v * (v^2-(3.0)*t^2) + d * u * (u^2-(3.0)*w^2) - b * w * (w^2-(3.0)*u^2) - 9.48) = 0;


solve;
	display a;
	display b;
	display c;
	display d;
	display t;
	display u;
	display v;
	display w;
display obj;
