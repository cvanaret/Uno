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
#   an exercize for L. Watson course on LANCELOT in the Spring 1993.

#   SIF input: T. Kuan, Virginia Tech., Spring 1993.

#   classification QLR2-AN-12-7

param x{1..6};
param a{1..6};
param c{1..7};

var t{1..6} >= 0;
var u{1..6} >= 0;

minimize f:
	sum {i in 1..6} (x[i]*t[i])^2;
subject to cons1:
	a[1]*t[1]-u[1]-c[1] = 0;
subject to cons2:
	a[2]*t[2]-u[2]-t[2]-u[2]-t[3]-u[3]-c[2] = 0;
subject to cons3:
	a[3]*t[3]-u[3]-t[3]-u[3]-t[1]-u[1]-t[6]-u[6]-t[4]-u[4]-t[5]-u[5]-t[2]-u[2]-c[3]=0;
subject to cons4:
	a[4]*t[4]-u[4]-t[4]-u[4]-t[1]-u[1]-t[5]-u[5]-t[6]-u[6]-c[4]=0;
subject to cons5:
	a[5]*t[5]-u[5]-t[1]-u[1]-c[5] = 0;
subject to cons6:
	a[6]*t[6] - u[6] + sum {i in 1..5} (-t[i]-u[i]) - c[6] = 0;
subject to cons7:
	sum {i in 1..6} (t[i]+u[i]) - c[7] = 0;

data;
param:	x	a:=	
1	1.502	1.8	
2	1.126	3.2	
3	0.815	6.1	
4	1.268	3.2	
5	1.502	1.8	
6	0.740	7.4;	
param c:=
1	11.0
2	3.0
3	20.0
4	17.0
5	9.0
6	20.0
7	126.1;

solve; display f; 
display t,u;
