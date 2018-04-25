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
#   W. Pulleyblank,
#   private communication, 1991.

#   classification LQR2-MY-36-66

param nnodes:=12;
param epsil := 0.0001;
param a{i in 1..nnodes, j in 1..i-1} default 0.0;

var x{i in 1..nnodes} := i*0.03;
var y{i in 1..nnodes} := i*0.0055555;
var r{i in 1..nnodes} >= 1, := i;

minimize f:
	sum {i in 1..nnodes} r[i];
subject to cons1{i in 2..nnodes, j in 1..i-1: a[i,j] = 1}:
	(r[i]+r[j])^2-(x[i]-x[j])^2-(y[i]-y[j])^2 = 0;
subject to cons2{i in 2..nnodes, j in 1..i-1: a[i,j] = 0}:
	(r[i]+r[j])^2-(x[i]-x[j])^2-(y[i]-y[j])^2 + epsil <= 0;
subject to cons3:
	x[1] = 0;
subject to cons4:
	y[1] = 0;
subject to cons5:
	y[2] = 0;

data;
param a:= 
2	1	1
7	1	1
3	2	1
4	2	1
4	3	1
5	4	1
6	4	1
6	5	1
11	5	1
7	6	1
8	7	1
9	7	1
9	8	1
10	8	1
10	9	1
11	10	1
12	10	1
12	11	1;

solve;
display f;
display x,y,r;
