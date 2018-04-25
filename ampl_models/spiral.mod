#
#**************************
# SET UP THE INITIAL DATA *
#**************************
#   Problem :
#   *********
#   A nonlinear minmax problem.
#   Source:
#   E. Polak, J.E. Higgins and D. Mayne,
#   "A barrier function for minmax problems",
#   Mathematical Programming, vol.54(2), pp. 155-176, 1992.
#   SIF input: Ph. Toint, April 1992.
#   classification LOR2-AN-3-2
#   Solution

	var x1 := 1.41831;
	var x2 := -4.79462;
	var u := 1.0;

minimize obj:
	u;

subject to c1:
	0 <= -(x1-(sqrt((x1*x1)+(x2*x2)))*cos((sqrt((x1*x1)+(x2*x2))))) * 
	(x1-(sqrt((x1*x1)+(x2*x2)))*cos((sqrt((x1*x1)+(x2*x2))))) - 0.0050*x1 * x1 - 
	0.0050*x2 * x2 + u;
subject to c2:
	0 <= -(x2-(sqrt((x1*x1)+(x2*x2)))*sin((sqrt((x1*x1)+(x2*x2))))) * 
	(x2-(sqrt((x1*x1)+(x2*x2)))*sin((sqrt((x1*x1)+(x2*x2))))) - 0.0050*x1 * x1 - 
	0.0050*x2 * x2 + u;

solve;
	display x1;
	display x2;
	display u;
display obj;
