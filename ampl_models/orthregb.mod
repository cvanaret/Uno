#
#**************************
# SET UP THE INITIAL DATA *
#**************************
#   Problem :
#   *********
#   An orthogonal regression problem.
#   The problem is to fit (orthogonally) an ellipse to a set of 6 points
#   in the 3D space. These points are compatible with this constraint.
#   Source:
#   M. Gulliksson,
#   "Algorithms for nonlinear Least-squares with Applications to
#   Orthogonal Regression",
#   UMINF-178.90, University of Umea, Sweden, 1990.
#   SIF input: Ph. Toint, June 1990.
#              correction by Ph. Shott, Jan 1995.
#   classification QQR2-AN-27-6
#   Parameters for the generation of the data points
#   Constants
#   Computed parameters
#   Construct the data points
#   Parameters of the ellipse
#   Projections of the data points onto the ellipse
#   Solution
	param a := 9.0;
	param b := 6.0;
	param c := 7.0;
	param cx := 0.5;
	param cy := 0.5;
	param cz := 0.5;
	param npts := 1 + (1 + (1 + (1 + (1 + (1)))));
	param xz := 0.5;
	param yz := 0.5;
	param zz := 0.5;
	param xd1 := (0.5) + (9.0);
	param yd1 := (0.5) + (9.0);
	param zd1 := 0.5;
	param xd2 := (0.5) + (6.0);
	param yd2 := (0.5) + (-1.0 * (6.0));
	param zd2 := 0.5;
	param xd3 := (0.5) + (-1.0 * (9.0));
	param yd3 := (0.5) + (-1.0 * (9.0));
	param zd3 := 0.5;
	param xd4 := (0.5) + (-1.0 * (6.0));
	param yd4 := (0.5) + (6.0);
	param zd4 := 0.5;
	param xd5 := 0.5;
	param yd5 := 0.5;
	param zd5 := (0.5) + (7.0);
	param xd6 := 0.5;
	param yd6 := 0.5;
	param zd6 := (0.5) + (-1.0 * (7.0));

	var h11 := 1.0;
	var h12;
	var h13;
	var h22 := 1.0;
	var h23;
	var h33 := 1.0;
	var g1;
	var g2;
	var g3;
	var x1 := 9.5;
	var y1 := 9.5;
	var z1 := 0.5;
	var x2 := 6.5;
	var y2 := -5.5;
	var z2 := 0.5;
	var x3 := -8.5;
	var y3 := -8.5;
	var z3 := 0.5;
	var x4 := -5.5;
	var y4 := 6.5;
	var z4 := 0.5;
	var x5 := 0.5;
	var y5 := 0.5;
	var z5 := 7.5;
	var x6 := 0.5;
	var y6 := 0.5;
	var z6 := -6.5;

minimize obj:
	(x1 - 9.5)*(x1 - 9.5) + (y1 - 9.5)*(y1 - 9.5) + (z1 - 0.5)*(z1 - 0.5) + (x2 - 
	6.5)*(x2 - 6.5) + (y2 + 5.5)*(y2 + 5.5) + (z2 - 0.5)*(z2 - 0.5) + (x3 + 
	8.5)*(x3 + 8.5) + (y3 + 8.5)*(y3 + 8.5) + (z3 - 0.5)*(z3 - 0.5) + (x4 + 
	5.5)*(x4 + 5.5) + (y4 - 6.5)*(y4 - 6.5) + (z4 - 0.5)*(z4 - 0.5) + (x5 - 
	0.5)*(x5 - 0.5) + (y5 - 0.5)*(y5 - 0.5) + (z5 - 7.5)*(z5 - 7.5) + (x6 - 
	0.5)*(x6 - 0.5) + (y6 - 0.5)*(y6 - 0.5) + (z6 + 6.5)*(z6 + 6.5);

subject to e1:
	h11 * x1 * x1 + 2.0*h12 * x1 * y1 + h22 * y1 * y1 - 2.0*g1 * x1 - 2.0*g2 * y1 + 
	2.0*h13 * x1 * z1 + 2.0*h23 * y1 * z1 + h33 * z1 * z1 - 2.0*g3 * z1 - 1.0 = 0;
subject to e2:
	h11 * x2 * x2 + 2.0*h12 * x2 * y2 + h22 * y2 * y2 - 2.0*g1 * x2 - 2.0*g2 * y2 + 
	2.0*h13 * x2 * z2 + 2.0*h23 * y2 * z2 + h33 * z2 * z2 - 2.0*g3 * z2 - 1.0 = 0;
subject to e3:
	h11 * x3 * x3 + 2.0*h12 * x3 * y3 + h22 * y3 * y3 - 2.0*g1 * x3 - 2.0*g2 * y3 + 
	2.0*h13 * x3 * z3 + 2.0*h23 * y3 * z3 + h33 * z3 * z3 - 2.0*g3 * z3 - 1.0 = 0;
subject to e4:
	h11 * x4 * x4 + 2.0*h12 * x4 * y4 + h22 * y4 * y4 - 2.0*g1 * x4 - 2.0*g2 * y4 + 
	2.0*h13 * x4 * z4 + 2.0*h23 * y4 * z4 + h33 * z4 * z4 - 2.0*g3 * z4 - 1.0 = 0;
subject to e5:
	h11 * x5 * x5 + 2.0*h12 * x5 * y5 + h22 * y5 * y5 - 2.0*g1 * x5 - 2.0*g2 * y5 + 
	2.0*h13 * x5 * z5 + 2.0*h23 * y5 * z5 + h33 * z5 * z5 - 2.0*g3 * z5 - 1.0 = 0;
subject to e6:
	h11 * x6 * x6 + 2.0*h12 * x6 * y6 + h22 * y6 * y6 - 2.0*g1 * x6 - 2.0*g2 * y6 + 
	2.0*h13 * x6 * z6 + 2.0*h23 * y6 * z6 + h33 * z6 * z6 - 2.0*g3 * z6 - 1.0 = 0;

solve;
	display h11;
	display h12;
	display h13;
	display h22;
	display h23;
	display h33;
	display g1;
	display g2;
	display g3;
	display x1;
	display y1;
	display z1;
	display x2;
	display y2;
	display z2;
	display x3;
	display y3;
	display z3;
	display x4;
	display y4;
	display z4;
	display x5;
	display y5;
	display z5;
	display x6;
	display y6;
	display z6;
display obj;
