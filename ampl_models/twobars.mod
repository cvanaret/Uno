#
#**************************
# SET UP THE INITIAL DATA *
#**************************
#   Problem:
#   ********
#   Structureal analysis of the simplest two bar scheme.  The structure has
#   the following simple symmetric shape
#                                *
#                               / \
#                              /   \
#                             /     \
#                           """     """
#   and a force is applied at the top node.  The unknown are the distance
#   of the left and right feet wrt to the projection of the top node and the
#   weight of the bars.
#   Source:
#   an example in a talk by W.K. Zhang and C. Fleury, LLN, 1994.
#   SIF input: Ph. Toint, November 1994
#   classification OOR2-MN-2-2
#   Solution

	var x1 >= 0.2 ,  <= 4.0 ,  := 1.0;
	var x2 >= 0.1 ,  <= 1.6 ,  := 1.0;

minimize obj:
	x1 * (sqrt((1.0+x2*x2)));

subject to cons1:
	0 >= 0.124*(sqrt((1.0+x2*x2))) * ((8.0/x1)+1.0/(x1*x2)) - 1.0;
subject to cons2:
	0 >= 0.124*(sqrt((1.0+x2*x2))) * ((8.0/x1)-1.0/(x1*x2)) - 1.0;

solve;
	display x1;
	display x2;
display obj;
