# AMPL Model by Hande Y. Benson
#
# Copyright (C) 2001 Princeton University
# All Rights Reserved
#
# Permission to use, copy, modify, and distribute this software and
# its documentation for any purpose and without fee is hereby
# granted, provided that the above copyright notice appear in all
# copies and that the copyright notice and this
# permission notice appear in all supporting documentation.                     

#   classification LQI2-RN-65-59

var Flow{1..24};
var Prod{1..6};
var Supply{1..6};
var Demand{1..9};
var Pi{1..20};

param Region{1..9};

param loflows{1..24};  
param hiflows;  
param hisupply{1..6};  
param loprods{1..6};  
param hiprods{1..6};  
param lopi{1..20}; 
param hipi{1..20};

minimize cost:
	2.28*Prod[1] + 2.28*Prod[2] + 2.28*Prod[3] + 1.68*Prod[4] + 
	1.68*Prod[5] + 1.68*Prod[6];
subject to node1:
	Flow[1] + Flow[2] - Supply[1] = 0;
subject to node2:
	-1*Flow[1] - Flow[2] + Flow[3] + Flow[4] - Supply[2] = 0;
subject to node3:
	-1*Flow[3] - Flow[4] + Flow[5] + Demand[1] = 0;
subject to node4:
	-1*Flow[5] - Flow[8] + Flow[9] = 0;
subject to node5:
	Flow[6] - Supply[3] = 0;
subject to node6:
	-1*Flow[6] + Flow[7] + Demand[2] = 0;
subject to node7:
	-1*Flow[7] + Flow[8] + Demand[3] = 0;
subject to node8:
	Flow[10] + Flow[11] - Supply[4] = 0;
subject to node9:
	-1*Flow[10] - Flow[11] + Flow[12] + Flow[13] = 0;
subject to node10:
	-1*Flow[12] - Flow[13] + Flow[14] + Flow[15] + Demand[4] = 0;
subject to node11:
	-1*Flow[14] - Flow[15] + Flow[16] + Flow[21] = 0;
subject to node12:
	-1*Flow[16] + Flow[17] + Demand[5] = 0;
subject to node13:
	-1*Flow[17] + Flow[18] - Supply[5] = 0;
subject to node14:
	-1*Flow[9] - Flow[18] + Flow[19] - Supply[6] = 0;
subject to node15:
	-1*Flow[19] + Flow[20] + Demand[6] = 0;
subject to node16:
	-1*Flow[20] + Demand[7] = 0;
subject to node17:
	-1*Flow[21] + Flow[22] =0;
subject to node18:
	-1*Flow[22] + Flow[23] = 0;
subject to node19:
	-1*Flow[23] + Flow[24] + Demand[8] = 0;
subject to node20:
	-1*Flow[24] + Demand[9] = 0;
subject to region1:
	-1*Demand[1] <= Region[1];
subject to region2:
	-1*Demand[2] <= Region[2];
subject to region3:
	-1*Demand[3] <= Region[3];
subject to region4:
	-1*Demand[4] <= Region[4];
subject to region5:
	-1*Demand[5] <= Region[5];
subject to region6:
	-1*Demand[6] <= Region[6];
subject to region7:
	-1*Demand[7] <= Region[7];
subject to region8:
	-1*Demand[8] <= Region[8];
subject to region9:
	-1*Demand[9] <= Region[9];
subject to prod1:
	Supply[1] <= Prod[1];
subject to prod2:
	Supply[2] <= Prod[2];
subject to prod3:
	Supply[3] <= Prod[3];
subject to prod4:
	Supply[4] <= Prod[4];
subject to prod5:
	Supply[5] <= Prod[5];
subject to prod6:
	Supply[6] <= Prod[6];
subject to arc1:
	(Flow[1]^2) -9.07027*(Pi[1] - Pi[2]) = 0;
subject to arc2:
	(Flow[2]^2) -9.07027*(Pi[1] - Pi[2]) = 0;
subject to arc3:
	(Flow[3]^2) -6.04685*(Pi[2] - Pi[3])  = 0;
subject to arc4:
	(Flow[4]^2) -6.04685*(Pi[2] - Pi[3]) = 0;
subject to arc5:
	(Flow[5]^2) -1.39543*(Pi[3] - Pi[4]) = 0;
subject to arc6:
	(Flow[6]^2) -0.100256*(Pi[5] - Pi[6]) = 0;
subject to arc7:
	(Flow[7]^2) -0.148655*(Pi[6] - Pi[7]) = 0;
subject to arc8:
	(Flow[8]^2) +0.226895*(Pi[4] - Pi[7]) = 0;
subject to arc9:
	(Flow[9]^2) -0.659656*(Pi[4] - Pi[14]) = 0;
subject to arc10:
	(Flow[10]^2) -7.25622*(Pi[8] - Pi[9]) >= 0;
subject to arc11:
	(Flow[11]^2) -0.10803*(Pi[8] - Pi[9]) >= 0;
subject to arc12:
	(Flow[12]^2) -1.81405*(Pi[9] - Pi[10]) = 0;
subject to arc13:
	(Flow[13]^2) -0.0270084*(Pi[9] - Pi[10]) = 0;
subject to arc14:
	(Flow[14]^2) -1.45124*(Pi[10] - Pi[11]) = 0;
subject to arc15:
	(Flow[15]^2) -0.0216067*(Pi[10] - Pi[11]) = 0;
subject to arc16:
	(Flow[16]^2) -0.863836*(Pi[11] - Pi[12]) = 0;
subject to arc17:
	(Flow[17]^2) -0.907027*(Pi[12] - Pi[13]) = 0;
subject to arc18:
	(Flow[18]^2) -7.25622*(Pi[13] - Pi[14]) = 0;
subject to arc19:
	(Flow[19]^2) -3.62811*(Pi[14] - Pi[15]) = 0;
subject to arc20:
	(Flow[20]^2) -1.45124*(Pi[15] - Pi[16]) = 0;
subject to arc21:
	(Flow[21]^2) -0.0514445*(Pi[11] - Pi[17]) = 0;
subject to arc22:
	(Flow[22]^2) -0.00641977*(Pi[17] - Pi[18]) >= 0;
subject to arc23:
	(Flow[23]^2) -0.00170320*(Pi[18] - Pi[19]) = 0;
subject to arc24:
	(Flow[24]^2) -0.0278190*(Pi[19] - Pi[20]) = 0;
subject to f{i in 1..24}: 
	loflows[i] <= Flow[i] <= hiflows;  
subject to s{i in 1..6}: 
	Supply[i] <= hisupply[i];  
subject to pr{i in 1..6}: 
	loprods[i] <= Prod[i] <= hiprods[i];  
subject to pi{i in 1..20}: 
	lopi[i] <= Pi[i] <= hipi[i];

data;
param Region:=
1	-3.91800
2	-4.03400
3	-5.25600
4	-6.36500
5	-2.12000
6	-6.84800
7	-15.6160
8	-0.222000
9	-1.91900;

param loflows:= 1 -220.12 2 -220.12 3 -220.12 4 -220.12 5 -220.12 6
-220.12 7 -220.12 8 -220.12 9 -220.12 10 0 11 0 12 -220.12 13 -220.12 14
-220.12 15 -220.12 16 -220.12 17 -220.12 18 -220.12 19 -220.12 20 -220.12
21 -220.12 22 0 23 -220.12 24 -220.12;  
param hiflows := 220.12;  
param hisupply:= 1 11.594 2 8.4 3 4.8 4 22.012 5 1.2 6 0.96; 
param loprods:= 1 8.87 2 0 3 0 4 20.344 5 0 6 0;  
param hiprods:= 1 11.594 2 8.4 3 4.8 4 22.012 5 1.2 6 0.96;  
param lopi:= 1 0 2 0 3 900 4 0 5 0 6 900 7 900 8 2500 9 0 10 900 11 0 12 
900 13 0 14 0 15 900 16 2500 17 0 18 0 19 625 20 625;  
param hipi:= 1 5929 2 5929 3 6400 4 6400 5 5929 6 6400 7 6400 8
4382.44 9 4382.44 10 4382.44 11 4382.44 12 4382.44 13 4382.44 14 4382.44
15 4382.44 16 4382.44 17 4382.44 18 3969 19 4382.44 20 4382.44;  
var Flow:= 1 5.797 2 5.797 3 9.997 4 9.997 5 16.076 6 4.8 7 0.766 8 -4.49 9
11.586 10 17.2404 11 2.10363 12 17.2404 13 2.10363 14 11.5676 15 1.41145
16 10.838 17 8.718 18 9.918 19 22.464 20 15.616 21 2.141 22 2.141 23 2.141
24 1.919;  

solve; 

