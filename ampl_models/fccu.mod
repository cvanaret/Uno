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

#   Source:
#   W. J. Korchinski, Profimatics, Inc,
#   325 Rolling Oaks Drive, Thousand Oaks, California, USA 91361-1200
#   Telephone: 1-805 496 6661, Fax: 1-805 373 5108

#   SIF input: W. Korchinski, Spring 1993.

#   classification SLR2-MN-19-8

param w{1..19};
param m{1..19};

var Feed := 1;
var Effluent := 1; 
var MF_ohd := 1; 
var HCN := 1;
var LCO := 1;
var HCO := 1;
var MF_btms := 1;
var Decant := 1;
var Dec_recy := 1;
var Off_gas := 1;
var DC4_feed := 1;
var DC3_feed := 1;
var DC4_btms := 1;
var Lean_oil := 1;
var Propane := 1;
var Butane := 1;
var C8spl_fd := 1;
var LCN := 1;
var MCN := 1;

minimize f:
    (Feed-m[1])^2/w[1]
    +(Effluent-m[2])^2/w[2]
    +(MF_ohd-m[3])^2/w[3]
    +(HCN-m[4])^2/w[4]
    +(LCO-m[5])^2/w[5]
    +(HCO-m[6])^2/w[6]
    +(MF_btms-m[7])^2/w[7]
    +(Decant-m[8])^2/w[8]
    +(Dec_recy-m[9])^2/w[9]
    +(Off_gas-m[10])^2/w[10]
    +(DC4_feed-m[11])^2/w[11]
    +(DC3_feed-m[12])^2/w[12]
    +(DC4_btms-m[13])^2/w[13]
    +(Lean_oil-m[14])^2/w[14]
    +(Propane-m[15])^2/w[15]
    +(Butane-m[16])^2/w[16]
    +(C8spl_fd-m[17])^2/w[17]
    +(LCN-m[18])^2/w[18]
    +(MCN-m[19])^2/w[19]
;
subject to cons1:
	Feed + Dec_recy - Effluent = 0;
subject to cons2:
	Effluent - MF_ohd - HCN - LCO - HCO - MF_btms = 0;
subject to cons3:
	MF_btms - Decant - Dec_recy = 0;
subject to cons4:
	MF_ohd + Lean_oil - Off_gas - DC4_feed = 0;
subject to cons5:
	DC4_feed - DC3_feed - DC4_btms = 0;
subject to cons6:
	DC4_btms - Lean_oil - C8spl_fd = 0;
subject to cons7:
	DC3_feed - Propane - Butane = 0;
subject to cons8:
	C8spl_fd - LCN - MCN = 0;

data; 
param w:= 
1 0.2 
2 1 
3 1 
4 0.33333333 
5 0.33333333
6 0.33333333 
7 1 
8 1 
9 1 
10 1 
11 1 
12 1 
13 1 
14 1 
15 0.33333333 
16 0.33333333 
17 1 
18 0.33333333 
19 0.33333333;

param m:= 
1 31 
2 36 
3 20 
4 3 
5 5 
6 3.5 
7 4.2 
8 0.9 
9 3.9
10 2.2 
11 22.8 
12 6.8 
13 19 
14 8.5 
15 2.2 
16 2.5 
17 10.8 
18 6.5 
19 6.5;

solve;
display f;
display Feed;
display Effluent;
display MF_ohd;
display HCN;
display LCO;
display HCO;
display MF_btms;
display Decant;
display Dec_recy;
display Off_gas;
display DC4_feed;
display DC3_feed;
display DC4_btms;
display Lean_oil;
display Propane;
display Butane;
display C8spl_fd;
display LCN;
display MCN;
