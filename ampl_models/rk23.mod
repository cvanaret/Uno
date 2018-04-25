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
#   Similar ideas for 4th and 5th order pairs are discussed in:
#   Hairer, Norsett and Wanner, Solving Ordinary Differential
#   Equations I, Springer 1980, page 158 ff.

#   SIF input: S. Leyffer, January 1997.

#   classification LOR2-RN-17-11

var C2 := 1;
var A21 := 1;
var C3 := 0.5;
var A31 := 0.25;
var A32 := 0.25;
var B1 := 0.5;
var B2 := 0.5;
var B3 := 0;
var BB1 := 1/6;
var BB2 := 1/6;
var BB3 := 4/6;
var TP1 >= 0;
var TM1 >= 0;
var TP2 >= 0;
var TM2 >= 0;
var TP3 >= 0;
var TM3 >= 0;

minimize f:
	TP1+TM1+TP2+TM2+TP3+TM3;
subject to cons1:
	A21-C2 = 0;
subject to cons2:
	A31+A32-C3 = 0;
subject to cons3:
	B1+B2+B3-1 = 0;
subject to cons4:
	BB1+BB2+BB3-1 = 0;
subject to cons5:
	B2*C2+B3*C3-0.5 = 0;
subject to cons6:
	BB2*C2+BB3*C3-0.5 = 0;
subject to cons7:
	BB2*C2^2+BB3*C3^2-1/3 = 0;
subject to cons8:
	BB3*A32*C2-1/6 = 0;
subject to cons9:
	TP1-TM2-1+4*BB2*C2^3+4*BB3*C3^3 = 0;
subject to cons10:
	TP2-TM2-1+8*BB3*C3*A32*C2 = 0;
subject to cons11:
	TP3-TM3-1+12*BB3*A32*C2^2 = 0;

solve;
display f;
display C2; 
display A21;
display C3;
display A31;
display A32;
display B1;
display B2;
display B3;
display BB1;
display BB2;
display BB3;
display TP1;
display TM1;
display TP2;
display TM2;
display TP3;
display TM3;

