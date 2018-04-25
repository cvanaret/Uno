#
#**************************
# SET UP THE INITIAL DATA *
#**************************
#   Problem :
#   *********
#   Source: a very simple box-constrained quadratic
#   SIF input: Nick Gould, June 1996
#   classification QUR2-AN-2-0
#   Solution
	param one := 1.0;
	param theta := 1.0e10;
	param beta := (((1.0) / (1.0e10)) * ((1.0) / (1.0e10))) * ((1.0) / (1.0e10));
	param beta2 := ((1.0) / (1.0e10)) * ((1.0) / (1.0e10));

	var x1 := 1.0e-30;
	var x2 := 1.0;

minimize obj:
	0.5*(1.0e10*x1)*(1.0e10*x1) + 0.5*(x2)*(x2);


solve;
	display x1;
	display x2;
display obj;
