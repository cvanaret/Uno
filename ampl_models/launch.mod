#
#**************************
# SET UP THE INITIAL DATA *
#**************************
#   Problem:
#   ********
#   The objective function to be minimized represents the total cost of
#   the development and launching of a 3 stages space launching vehicle.
#   Constraints are imposed on physical interrelations between the variables
#   and performance.
#   The problem is highly non-convex. 
#   Source:
#   B. Rush, J. Bracken and G. McCormick,
#   "A nonliner programming model for launch vehicle design and costing",
#   Operations Research, pp. 185-210, 1967.
#   SIF input: P. Driscoll, Virginia Tech., April 1993,
#              corrected and simplified by Ph. L. Toint, May 1993.
#   classification OOR2-MY-25-28
#   Cost of stage 1
#   Cost of stage 2
#   Cost of stage 3
#   Instrumentation cost
#   Launch operating costs
#XN LAUN      'SCALE'   0.039215686
#   Relations between airframe weights and inert weights
#   Definition of stage mass fractions
#   Relation between stage thrust and single engine thrust
#   Constraints on structural integrity (stage weight vs propellant weight)
#   Constraints expression the ratio of thrust to initial weight for a given 
#   payload
#   Constraints on the stage mass fraction of the 3 stages
#   Constraints on the specific impulse
#   Total vehicle launch constraint
#   Objective function
#   Constraints
#   The starting point is very close to feasible
#   Elements for stage 1 modelling
#   a) airframe R&D
#   b) LOX/RP propulsion R&D
#   c) airframe production unit
#   d) LOX/RP propulsion production of 5 engines
#   Elements for stage 2 modelling
#   a) airframe R&D
#   b) LOX/RP propulsion R&D
#   c) airframe production unit
#   d) LOX/RP propulsion production of 5 engines
#   Elements for stage 3 modelling
#   a) airframe R&D
#   b) LOX/RP propulsion R&D
#   c) airframe production unit
#   d) LOX/RP propulsion production of 1 engine
#   Elements for constraints on the stage mass fraction
#   Constraints on specific impulse

var aw1 >= 1.0e-8 ,  <= 10000.0 ,  := 68.0;
var iw1 >= 1.0e-8 ,  <= 10000.0 ,  := 136.0;
var mf1 >= 0.25 ,  <= 0.3 ,  := 0.29988744;
var tt1 >= 1.0e-8 ,  <= 10000.0 ,  := 3733.0;
var pw1 >= 1.0e-8 ,  <= 10000.0 ,  := 2177.0;
var et1 >= 1.0e-8 ,  <= 10000.0 ,  := 746.6;
var s1l >= 125.0 ,  <= 150.0 ,  := 125.0;
var aw2 >= 1.0e-8 ,  <= 10000.0 ,  := 28.2;
var iw2 >= 1.0e-8 ,  <= 10000.0 ,  := 47.0;
var mf2 >= 0.24 ,  <= 0.29 ,  := 0.28939109;
var tt2 >= 1.0e-8 ,  <= 10000.0 ,  := 480.0;
var pw2 >= 1.0e-8 ,  <= 10000.0 ,  := 566.0;
var et2 >= 1.0e-8 ,  <= 10000.0 ,  := 96.0;
var s2l >= 75.0 ,  <= 100.0 ,  := 75.0;
var aw3 >= 1.0e-8 ,  <= 10000.0 ,  := 11.2;
var iw3 >= 1.0e-8 ,  <= 10000.0 ,  := 16.0;
var mf3 >= 0.16 ,  <= 0.21 ,  := 0.20980926;
var tt3 >= 1.0e-8 ,  <= 10000.0 ,  := 129.0;
var pw3 >= 1.0e-8 ,  <= 10000.0 ,  := 145.0;
var et3 >= 1.0e-8 ,  <= 10000.0 ,  := 129.0;
var s3l >= 50.0 ,  <= 70.0 ,  := 50.0;
var inw >= 2.5 ,  <= 4.0 ,  := 2.5;
var bt1 >= 1.0e-8 ,  <= 10000.0 ,  := 155.0;
var bt2 >= 1.0e-8 ,  <= 10000.0 ,  := 314.0;
var bt3 >= 1.0e-8 ,  <= 10000.0 ,  := 403.0;

minimize obj:
	((5272.77*(aw1^1.2781) * (iw1^-0.1959) * (mf1^2.4242) * (tt1^0.38745) * (pw1^0.9904) + 160.909*(et1 / 1000.0 ) ^-0.146 + 282.874*(et1 / 1000.0 ) ^0.648 + 0.64570846*(aw1^0.3322) * (mf1^-1.5935) * (pw1^0.2363) * (s1l^0.1079) + 31.136196*(et1 / 1000.0 ) ^0.736 + 12.092112*(et1 / 1000.0 ) ^-0.229 + 31.136196*(et2 / 1000.0 ) ^0.736 + 12.092112*(et2 / 1000.0 ) ^-0.229 + 2.5870000000000005e-4*et1 - 247.963)/1.0e8) + ((5272.77*(aw2^1.2781) * (iw2^-0.1959) * (mf2^2.4242) * (tt2^0.38745) * (pw2^0.9904) + 160.909*(et2 / 1000.0 ) ^-0.146 + 282.874*(et2 / 1000.0 ) ^0.648 + 0.64570846*(aw2^0.3322) * (mf2^-1.5935) * (pw2^0.2363) * (s2l^0.1079) + 2.5870000000000005e-4*et2 - 247.963)/1.0e8) + ((5272.77*(aw3^1.2781) * (iw3^-0.1959) * (mf3^2.4242) * (tt3^0.38745) * (pw3^0.9904) + 181.806*(et3 / 1000.0 ) ^0.539 + 232.57*(et3 / 1000.0 ) ^0.772 + 0.49783215*(aw3^0.3322) * (mf3^-1.5935) * (pw3^0.2363) * (s3l^0.1079) - 0.22424514*(et3 / 100.0 ) ^-1.33 + 20.708238*(et3 / 100.0 ) ^0.498 + 0.001958*et3 - 32.591)/1.0e8) + ((47.040096*inw - 35.5)/1.0e8) + (((0.0030*pw1 + 0.0030*pw2 + 0.0030*pw3)^0.460)/3.9215686e7);

subject to obj_bnd:
	0.0 <= ((5272.77*(aw1^1.2781) * (iw1^-0.1959) * (mf1^2.4242) * (tt1^0.38745) * (pw1^0.9904) + 160.909*(et1 / 1000.0 ) ^-0.146 + 282.874*(et1 / 1000.0 ) ^0.648 + 0.64570846*(aw1^0.3322) * (mf1^-1.5935) * (pw1^0.2363) * (s1l^0.1079) + 31.136196*(et1 / 1000.0 ) ^0.736 + 12.092112*(et1 / 1000.0 ) ^-0.229 + 31.136196*(et2 / 1000.0 ) ^0.736 + 12.092112*(et2 / 1000.0 ) ^-0.229 + 2.5870000000000005e-4*et1 - 247.963)/1.0e8) + ((5272.77*(aw2^1.2781) * (iw2^-0.1959) * (mf2^2.4242) * (tt2^0.38745) * (pw2^0.9904) + 160.909*(et2 / 1000.0 ) ^-0.146 + 282.874*(et2 / 1000.0 ) ^0.648 + 0.64570846*(aw2^0.3322) * (mf2^-1.5935) * (pw2^0.2363) * (s2l^0.1079) + 2.5870000000000005e-4*et2 - 247.963)/1.0e8) + ((5272.77*(aw3^1.2781) * (iw3^-0.1959) * (mf3^2.4242) * (tt3^0.38745) * (pw3^0.9904) + 181.806*(et3 / 1000.0 ) ^0.539 + 232.57*(et3 / 1000.0 ) ^0.772 + 0.49783215*(aw3^0.3322) * (mf3^-1.5935) * (pw3^0.2363) * (s3l^0.1079) - 0.22424514*(et3 / 100.0 ) ^-1.33 + 20.708238*(et3 / 100.0 ) ^0.498 + 0.001958*et3 - 32.591)/1.0e8) + ((47.040096*inw - 35.5)/1.0e8) + (((0.0030*pw1 + 0.0030*pw2 + 0.0030*pw3)^0.460)/3.9215686e7);
subject to sgth1:
	2.0*aw1 - iw1 = 0;
subject to sgth3:
	0.6*iw2 - aw2 = 0;
subject to sgth5:
	0.7*iw3 - aw3 = 0;
subject to sgth2:
	5.0*et1 - tt1 = 0;
subject to sgth4:
	5.0*et2 - tt2 = 0;
subject to sgth6:
	tt3 - et3 = 0;
subject to sgsi1a:
	0 <= pw1 - 12.0*iw1;
subject to sgsi1b:
	0 >= pw1 - 17.0*iw1;
subject to sgsi2a:
	0 <= pw2 - 10.0*iw2;
subject to sgsi2b:
	0 >= pw2 - 13.0*iw2;
subject to sgsi3a:
	0 <= pw3 - 7.0*iw3;
subject to sgsi3b:
	0 >= pw3 - 10.0*iw3;
subject to ttiw1a:
	0 <= tt1 - 1.2*iw1 - 1.2*pw1 - 1.2*iw2 - 1.2*pw2 - 1.2*iw3 - 1.2*pw3 - 1.2*inw - 24.0;
subject to ttiw1b:
	0 >= tt1 - 1.4*iw1 - 1.4*pw1 - 1.4*iw2 - 1.4*pw2 - 1.4*iw3 - 1.4*pw3 - 1.4*inw - 28.0;
subject to ttiw2a:
	0 <= tt2 - 0.6*iw2 - 0.6*pw2 - 0.6*iw3 - 0.6*pw3 - 0.6*inw - 12.0;
subject to ttiw2b:
	0 >= tt2 - 0.75*iw2 - 0.75*pw2 - 0.75*iw3 - 0.75*pw3 - 0.75*inw - 15.0;
subject to ttiw3a:
	0 <= tt3 - 0.7*iw3 - 0.7*pw3 - 0.7*inw - 14.0;
subject to ttiw3b:
	0 >= tt3 - 0.9*iw3 - 0.9*pw3 - 0.9*inw - 18.0;
subject to smf1:
	(mf1) * (inw) + 20.0*mf1 - iw1 - iw2 - pw2 - iw3 - pw3 - inw - 20.0 = 0;
subject to smf2:
	(mf2) * (inw) + 20.0*mf2 - iw2 - iw3 - pw3 - inw - 20.0 = 0;
subject to smf3:
	(mf3) * (inw) + 20.0*mf3 - iw3 - inw - 20.0 = 0;
subject to si1a:
	0 <= tt1 * bt1 - 240.0*pw1;
subject to si1b:
	0 >= tt1 * bt1 - 290.0*pw1;
subject to si2a:
	0 <= tt2 * bt2 - 240.0*pw2;
subject to si2b:
	0 >= tt2 * bt2 - 290.0*pw2;
subject to si3a:
	0 <= tt3 * bt3 - 340.0*pw3;
subject to si3b:
	0 >= tt3 * bt3 - 375.0*pw3;
subject to glgcon:
	0 <= -32.0*(tt1 * bt1 * (log(mf1)) ) /pw1 - 32.0*(tt2 * bt2 * (log(mf2)) ) /pw2 - 32.0*(tt3 * bt3 * (log(mf3)) ) /pw3 - 35000.0 <= 15000.0;

solve;
display aw1;
display iw1;
display mf1;
display tt1;
display pw1;
display et1;
display s1l;
display aw2;
display iw2;
display mf2;
display tt2;
display pw2;
display et2;
display s2l;
display aw3;
display iw3;
display mf3;
display tt3;
display pw3;
display et3;
display s3l;
display inw;
display bt1;
display bt2;
display bt3;
display obj;
