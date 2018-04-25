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
#   J.P. Rasson, Private communication, 1996.

#   SIF input: Ph. Toint, Nov 1996.

#   classification SBR2-MN-61-0

param LGSG := 11;
param LGTR := 40;

param PIC := 3.0000000000;

param TR{1..LGTR};
param SSG{1..LGSG};
param CC{1..40} := 0.0;

var c{i in 1..LGTR} := CC[i], >= 0.0;
var sg{i in 1..LGSG} := SSG[i], >= 0.0, <= PIC;
var x{k in 1..LGTR, i in 1..LGSG} = if k-i+1 <= 0 then 0 else
sg[i]*c[k-i+1];
minimize f:
	sum {k in 1..LGTR} (sum {i in 1..LGSG} x[k,i] 
	-TR[k])^2;

data;
param TR:=
1 0.0000000000 2 0.0000000000 3 1.600000E-03 4 5.400000E-03 5
7.020000E-02 6 0.1876000000 7 0.3320000000 8 0.7640000000 9
0.9320000000 10 0.8120000000 11 0.3464000000 12 0.2064000000 13
8.300000E-02 14 3.400000E-02 15 6.179999E-02 16 1.2000000000 17
1.8000000000 18 2.4000000000 19 9.0000000000 20 2.4000000000 21
1.8010000000 22 1.3250000000 23 7.620000E-02 24 0.2104000000 25
0.2680000000 26 0.5520000000 27 0.9960000000 28 0.3600000000 29
0.2400000000 30 0.1510000000 31 2.480000E-02 32 0.2432000000 33
0.3602000000 34 0.4800000000 35 1.8000000000 36 0.4800000000 37
0.3600000000 38 0.2640000000 39 6.000000E-03 40 6.000000E-03;

param SSG:=
1 1.000000E-02 2 2.000000E-02 3 0.4000000000 4 0.6000000000 5
0.8000000000 6 3.0000000000 7 0.8000000000 8 0.6000000000 9
0.4400000000 10 1.000000E-02 11 1.000000E-02;

solve;
display f;
display c,sg;

