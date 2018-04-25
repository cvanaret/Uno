#
#**************************
# SET UP THE INITIAL DATA *
#**************************
#   Problem:
#   ********
#   A nonlinear least-squares problem.  This problem arises in measuring
#   angles and distances to a vibrating beam using a laser-Doppler
#   velocimeter.
#   Source:
#   an exercize for L. Watson course on LANCELOT in the Spring 1993.
#   SIF input: B. E. Lindholm, Virginia Tech., Spring 1993.
#   classification SBR2-MN-3-0
	param zero := 0;
	param p := 16;
	param realp := 16.0;
	param y0 := 21.158931;
	param y1 := 17.591719;
	param y2 := 14.046854;
	param y3 := 10.519732;
	param y4 := 7.0058392;
	param y5 := 3.5007293;
	param y6 := 0.0;
	param y7 := -3.5007293;
	param y8 := -7.0058392;
	param y9 := -10.519732;
	param y10 := -14.046854;
	param y11 := -17.591719;
	param y12 := -21.158931;
	param y13 := -24.753206;
	param y14 := -28.379405;
	param y15 := -32.042552;
	param y16 := -35.747869;
	param index := 16;

	var alpha := 0.6;
	var beta := -0.6;
	var dist >= 0.0 ,  := 20.0;

minimize obj:
	(dist * (tan(alpha*(1.0-(0.0/16.0))+beta*(0.0/16.0))) - 21.158931)*(dist * 
	(tan(alpha*(1.0-(0.0/16.0))+beta*(0.0/16.0))) - 21.158931) + (dist * 
	(tan(alpha*(1.0-(1.0/16.0))+beta*(1.0/16.0))) - 17.591719)*(dist * 
	(tan(alpha*(1.0-(1.0/16.0))+beta*(1.0/16.0))) - 17.591719) + (dist * 
	(tan(alpha*(1.0-(2.0/16.0))+beta*(2.0/16.0))) - 14.046854)*(dist * 
	(tan(alpha*(1.0-(2.0/16.0))+beta*(2.0/16.0))) - 14.046854) + (dist * 
	(tan(alpha*(1.0-(3.0/16.0))+beta*(3.0/16.0))) - 10.519732)*(dist * 
	(tan(alpha*(1.0-(3.0/16.0))+beta*(3.0/16.0))) - 10.519732) + (dist * 
	(tan(alpha*(1.0-(4.0/16.0))+beta*(4.0/16.0))) - 7.0058392)*(dist * 
	(tan(alpha*(1.0-(4.0/16.0))+beta*(4.0/16.0))) - 7.0058392) + (dist * 
	(tan(alpha*(1.0-(5.0/16.0))+beta*(5.0/16.0))) - 3.5007293)*(dist * 
	(tan(alpha*(1.0-(5.0/16.0))+beta*(5.0/16.0))) - 3.5007293) + (dist * 
	(tan(alpha*(1.0-(6.0/16.0))+beta*(6.0/16.0))))*(dist * 
	(tan(alpha*(1.0-(6.0/16.0))+beta*(6.0/16.0)))) + (dist * 
	(tan(alpha*(1.0-(7.0/16.0))+beta*(7.0/16.0))) + 3.5007293)*(dist * 
	(tan(alpha*(1.0-(7.0/16.0))+beta*(7.0/16.0))) + 3.5007293) + (dist * 
	(tan(alpha*(1.0-(8.0/16.0))+beta*(8.0/16.0))) + 7.0058392)*(dist * 
	(tan(alpha*(1.0-(8.0/16.0))+beta*(8.0/16.0))) + 7.0058392) + (dist * 
	(tan(alpha*(1.0-(9.0/16.0))+beta*(9.0/16.0))) + 10.519732)*(dist * 
	(tan(alpha*(1.0-(9.0/16.0))+beta*(9.0/16.0))) + 10.519732) + (dist * 
	(tan(alpha*(1.0-(10.0/16.0))+beta*(10.0/16.0))) + 14.046854)*(dist * 
	(tan(alpha*(1.0-(10.0/16.0))+beta*(10.0/16.0))) + 14.046854) + (dist * 
	(tan(alpha*(1.0-(11.0/16.0))+beta*(11.0/16.0))) + 17.591719)*(dist * 
	(tan(alpha*(1.0-(11.0/16.0))+beta*(11.0/16.0))) + 17.591719) + (dist * 
	(tan(alpha*(1.0-(12.0/16.0))+beta*(12.0/16.0))) + 21.158931)*(dist * 
	(tan(alpha*(1.0-(12.0/16.0))+beta*(12.0/16.0))) + 21.158931) + (dist * 
	(tan(alpha*(1.0-(13.0/16.0))+beta*(13.0/16.0))) + 24.753206)*(dist * 
	(tan(alpha*(1.0-(13.0/16.0))+beta*(13.0/16.0))) + 24.753206) + (dist * 
	(tan(alpha*(1.0-(14.0/16.0))+beta*(14.0/16.0))) + 28.379405)*(dist * 
	(tan(alpha*(1.0-(14.0/16.0))+beta*(14.0/16.0))) + 28.379405) + (dist * 
	(tan(alpha*(1.0-(15.0/16.0))+beta*(15.0/16.0))) + 32.042552)*(dist * 
	(tan(alpha*(1.0-(15.0/16.0))+beta*(15.0/16.0))) + 32.042552) + (dist * 
	(tan(alpha*(1.0-(16.0/16.0))+beta*(16.0/16.0))) + 35.747869)*(dist * 
	(tan(alpha*(1.0-(16.0/16.0))+beta*(16.0/16.0))) + 35.747869);


solve;
	display alpha;
	display beta;
	display dist;
display obj;
