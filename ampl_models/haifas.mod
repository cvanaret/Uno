#
#**************************
# SET UP THE INITIAL DATA *
#**************************
#   Problem :
#   *********
#  Truss Topology Design (t6-9)
#  Source: M. Tsibulevsky, Optimization Laboratory,
#          Faculty of Industrial Engineering, Technion,
#          Haifa, 32000, Israel.
#  SIF input: Conn, Gould and Toint, May, 1992
#             minor correction by Ph. Shott, Jan 1995.
#   classification LQR2-AN-13-9
#  2 * Number of nodes
#  Number of potential bars
#   Define constants
	param n := 12;
	param m := 9;
	param j := 9;
	param i1 := 6;
	param i2 := 6;
	param l := 21;

	var z;
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

minimize obj:
	z;

subject to g1:
	0 >= 10.0*0.5 * x4 * x4 - z - x10;
subject to g2:
	0 >= 6.4*0.5 * x5 * x5 + 6.4*0.5 * x5 * x11 + 1.6*0.5 * x11 * x11 - z - x10;
subject to g3:
	0 >= 40.0*0.5 * x10 * x10 - 80.0*0.5 * x10 * x11 + 40.0*0.5 * x11 * x11 - z - 
	x10;
subject to g4:
	0 >= 6.4*0.5 * x4 * x4 - 6.4*0.5 * x4 * x10 + 1.6*0.5 * x10 * x10 - z - x10;
subject to g5:
	0 >= 10.0*0.5 * x5 * x5 - z - x10;
subject to g6:
	0 >= 6.4*0.5 * x6 * x6 + 6.4*0.5 * x6 * x12 + 1.6*0.5 * x12 * x12 - z - x10;
subject to g7:
	0 >= 40.0*0.5 * x11 * x11 - 80.0*0.5 * x11 * x12 + 40.0*0.5 * x12 * x12 - z - 
	x10;
subject to g8:
	0 >= 6.4*0.5 * x5 * x5 - 6.4*0.5 * x5 * x11 + 1.6*0.5 * x11 * x11 - z - x10;
subject to g9:
	0 >= 10.0*0.5 * x6 * x6 - z - x10;

solve;
	display z;
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
display obj;
