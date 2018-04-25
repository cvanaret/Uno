#
#**************************
# SET UP THE INITIAL DATA *
#**************************
#   Problem :
#   *********
#   Source:  problem 207 in
#   A.R. Buckley,
#   "Test functions for unconstrained minimization",
#   TR 1989CS-3, Mathematics, statistics and computing centre,
#   Dalhousie University, Halifax (CDN), 1989.
#   SIF input: Ph. Toint, Dec 1989.
#   classification NOR2-AN-2-2
#   Solution

	var x1;
	var x2;

minimize obj: 0;
subject to cons1:
	((x1-x2*x2) * (x1-sin(x2))) = 0;
subject to cons2:
	(((cos(x2))-x1) * (x2-cos(x1))) = 0;


solve;
	display x1;
	display x2;
display obj;
