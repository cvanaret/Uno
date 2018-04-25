#
#**************************
# SET UP THE INITIAL DATA *
#**************************
#   Problem :
#   *********
#   A small gas network problem.
#   SIF input: Sybille Schachler, Oxford, August 1992.
#              minor correction by Ph. Shott, Jan 1995.
	param n := 7;
	param m := 4;
	param demand := -1000.0;
	param pmax1 := 914.73;
	param pmax2 := 904.73;
	param k := -0.597053452;

	var p1 <= 914.73 ,  := 965;
	var p2 := 965;
	var p3 <= 904.73 ,  := 965;
	var p4 := 965;
	var p5 <= 904.73 ,  := 965;
	var p6 := 965;
	var p7 <= 914.73 ,  := 965;
	var q1 := 100.0;
	var f1 := 1000.0;
	var q2 := 100.0;
	var f2 := 1000.0;
	var q3 := -100.0;
	var f3 := 1000.0;
	var q4 := -100.0;
	var f4 <= 400.0 ,  := 1000.0;

minimize obj:
	 - p1 - p2 - p3 - p4 - p5 - p6 - p7;

subject to pan1:
	p1 * abs ( p1 )  - p2 * abs ( p2 )  - 0.597053452*q1 * abs ( q1 ) ^0.8539 = 0;
subject to pan2:
	p3 * abs ( p3 )  - p4 * abs ( p4 )  - 0.597053452*q2 * abs ( q2 ) ^0.8539 = 0;
subject to pan3:
	p4 * abs ( p4 )  - p5 * abs ( p5 )  - 0.597053452*q3 * abs ( q3 ) ^0.8539 = 0;
subject to pan4:
	p6 * abs ( p6 )  - p7 * abs ( p7 )  - 0.597053452*q4 * abs ( q4 ) ^0.8539 = 0;
subject to mbal1:
	q1 - f3 = 0;
subject to mbal2:
	-q1 + f1 = 0;
subject to mbal3:
	q2 - f1 = 0;
subject to mbal4:
	-q2 + q3 + 1000.0 = 0;
subject to mbal5:
	-q3 - f2 = 0;
subject to mbal6:
	q4 + f2 = 0;
subject to mbal7:
	-q4 - f4 = 0;

solve;
	display p1;
	display p2;
	display p3;
	display p4;
	display p5;
	display p6;
	display p7;
	display q1;
	display f1;
	display q2;
	display f2;
	display q3;
	display f3;
	display q4;
	display f4;
display obj;
