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

#   classification NOR2-MN-20-14

var L0 <= 100, >= 0, := 1.5000D-01;
var N <= 100, >= 0, := 2.4079D+01;
var F <= 30, >= 0, := 9.2459D-15;
var K <= 100, >= 0 , := 0.0000D+00;
var LB <= 50, >= 0, := 0.0000D+00;
var L <= 50, >= 0, := 1.5000D-01;
var DE <= 30, >= 0, := 6.8120D+00;
var DI <= 30, >= 0, := 6.6120D+00;
var TO <= 800, >= 0, :=0.0000D+00 ;
var TOB <= 800, >= 0, := 0.0000D+00;
var NU <= 50, >= 0.5, := 2.2079D+01;
var D <= 10, >= 0.1, := 1.0000D-01;
var P <= 20, >= 0, := 6.5268D-01;
var E <= 10, >= 0, := 5.5268D-01;
var P0 <= 1000, >= 1, := 6.5887D+02;
var G <= 80000, >= 40000, := 6.5887D+04;
var DM <= 30, >= 0.1, := 6.7120D+00;
var FR <= 50, >= 0, := 1.5000D-01;
var TOLIM <= 1000, >= 100, := 1.0000D+02;
var TOBLIM <= 1000, >= 100, := 1.0000D+02;

minimize f: 0;
subject to cons1:
	(-F)=0;
subject to cons2:
	(-K)=0;
subject to cons3:
	(DE-D-DM)=0;
subject to cons4:
	(DI+D-DM)=0;
subject to cons5:
	(D-P+E)=0;
subject to cons6:
	(NU-N+2)=0;
subject to cons7:
	(1.5*D-L0)=0;
subject to cons8:
	(L-LB-FR)=0;
subject to cons9:
	(LB)=0;
subject to cons10:
	(L-L0+F)=0;
subject to cons11:
	(TO)=0;
subject to cons12:
	(TOB) = 0;
subject to cons13:
	TO-TOLIM <= 0;
subject to cons14:
	TOB-TOBLIM <= 0;

solve;
display f;
display L0;
display N;
display F;
display K;
display LB;
display L;
display DE;
display DI;
display TO;
display TOB;
display NU;
display D;
display P;
display E;
display P0;
display G;
display DM;
display FR;
display TOLIM;
display TOBLIM;

