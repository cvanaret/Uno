#
#**************************
# SET UP THE INITIAL DATA *
#**************************
#   Problem :
#   *********
#   GROWTH problem in 3 variables
#   Fit the observed growth g(n) from Gaussian Elimination
#   with complete pivoting to a function of the form
#
#        U1 * n ** ( U2 + LOG(n) * U3 )
#   SIF input: Nick Gould, Nov, 1991.
#   classification NOR2-AN-3-12
#T  'DEFAULT' L2
#   Solution
	param n := 3;

	var u1 := 100.0;
	var u2;
	var u3;

minimize obj:
	(u1 * (8.0^(u2+(log(8.0))*u3)) - 8.0)^2 + (u1 * (9.0^(u2+(log(9.0))*u3)) - 
	8.4305)^2 + (u1 * (10.0^(u2+(log(10.0))*u3)) - 9.5294)^2 + (u1 * 
	(11.0^(u2+(log(11.0))*u3)) - 10.4627)^2 + (u1 * (12.0^(u2+(log(12.0))*u3)) - 
	12.0)^2 + (u1 * (13.0^(u2+(log(13.0))*u3)) - 13.0205)^2 + (u1 * 
	(14.0^(u2+(log(14.0))*u3)) - 14.5949)^2 + (u1 * (15.0^(u2+(log(15.0))*u3)) - 
	16.1078)^2 + (u1 * (16.0^(u2+(log(16.0))*u3)) - 18.0596)^2 + (u1 * 
	(18.0^(u2+(log(18.0))*u3)) - 20.4569)^2 + (u1 * (20.0^(u2+(log(20.0))*u3)) - 
	24.25)^2 + (u1 * (25.0^(u2+(log(25.0))*u3)) - 32.9863)^2;


solve;
	display u1;
	display u2;
	display u3;
display obj;
