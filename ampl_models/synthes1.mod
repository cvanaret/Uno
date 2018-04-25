#
#**************************
# SET UP THE INITIAL DATA *
#**************************
#   Problem :
#   *********
#   Source: Test problem 1 (Synthesis of processing system) in 
#   M. Duran & I.E. Grossmann,
#   "An outer approximation algorithm for a class of mixed integer nonlinear
#    programs", Mathematical Programming 36, pp. 307-339, 1986.
#   SIF input: S. Leyffer, October 1997
#   classification OOR2-AN-6-6

	var x1 >= 0.0 ,  <= 2.0;
	var x2 >= 0.0 ,  <= 2.0;
	var x3 >= 0.0 ,  <= 1.0;
	var y1 >= 0.0 ,  <= 1.0;
	var y2 >= 0.0 ,  <= 1.0;
	var y3 >= 0.0 ,  <= 1.0;

minimize obj:
	 - 18.0*log ( x2 + 1.0 )  - 19.2*log ( x1 - x2 + 1.0 )  + 5.0*y1 + 6.0*y2 + 
	8.0*y3 + 10.0*x1 - 7.0*x3 + 10.0;

subject to n1:
	0 <= 0.8*log ( x2 + 1.0 )  + 0.96*log ( x1 - x2 + 1.0 )  - 0.8*x3;
subject to n2:
	0 <= log ( x2 + 1.0 )  + 1.2*log ( x1 - x2 + 1.0 )  - x3 - 2.0*y3 + 2.0;
subject to l3:
	0 >= x2 - x1;
subject to l4:
	0 >= x2 - 2.0*y1;
subject to l5:
	0 >= -x2 + x1 - 2.0*y2;
subject to l6:
	0 >= y1 + y2 - 1.0;

solve;
	display x1;
	display x2;
	display x3;
	display y1;
	display y2;
	display y3;
display obj;
