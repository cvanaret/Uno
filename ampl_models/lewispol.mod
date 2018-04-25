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
#   A. Lewis, private communication.

#   SIF input: A.R. Conn and Ph. Toint, March 1990.

#   classification QOR2-AN-6-9

param N:=6;
param DEG := 3;
param PEN := 1D4;

param c{i in 0..DEG-1, j in i..N-1} := if (i == 0) then 1 else
				       c[i-1,j]*j;
param ct{i in 0..DEG-1} := if (i==0) then -1 else
			   ct[i-1]*(N-i+1);
param a_init{0..N-1};
var a{i in 0..N-1} >= -10, <= 10, := a_init[i];

minimize f:
	sum {j in 0..N-1} a[j]^2;
subject to cons1:
	sum {j in 0..N-1} a[j]*c[0,j] - ct[0] = 0;
subject to cons2{i in 1..DEG-1}:
	sum {j in i..N-1} a[j]*c[i,j] - ct[i] = 0;
subject to cons3{j in 0..N-1}:
	(a[j]^3-a[j]) / PEN = 0;

data;
param a_init:=
0	-1
1	1
2	1
3	0
4	1
5	-1;

solve;
display f;
display a;
