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

#   Source: Optimal Design of Multiproduct Batch Plant
#   G.R. Kocis & I.E. Grossmann,
#   "Global OPtimization of Nonconvex Mixed Integer Nonlinear Programmming
#    (MINLP) problems in Process Synthesis", Indust. Engng. Chem. Res.,
#   No. 27, pp 1407--1421, 1988.

#   SIF input: S. Leyffer, October 1997

#   classification OOR2-AN-46-73

# modified to a continuous problem

param M := 6;
param N := 5;
param NU := 4;
param VL := log(300);
param VU := log(3000);
param H := 6000;
param TLO{1..5};
param TUP{1..5};
param BLO{1..5};
param BUP{1..5};
param Q{1..5};

param svals{1..N,1..M};
param tvals{1..N,1..M};

param S{i in 1..N, j in 1..M} := log(svals[i,j]);
param T{i in 1..N, j in 1..M} := log(tvals[i,j]);

param alpha := 250;
param beta := 0.6;

var n{1..M} <= log(NU), >= 0.0;
var v{1..M} <= VU, >= VL;
var b{i in 1..N} <= BUP[i], >= BLO[i];
var tl{i in 1..N} <= TUP[i], >= TLO[i];
var y{1..NU,1..M} <= 1.0, >= 0.0;

minimize f:
	sum {j in 1..M} alpha*exp(n[j]+beta*v[j]);
subject to cons1{i in 1..N, j in 1..M}:
	v[j] - b[i] >= S[i,j];
subject to cons2{i in 1..N, j in 1..M}:
	n[j] + tl[i] >= T[i,j];
subject to cons3:
	sum {i in 1..N} Q[i]*exp(tl[i] - b[i]) <= H;
subject to cons4{j in 1..M}:
	sum {k in 1..NU} log(k)*y[k,j] - n[j] = 0;
subject to cons5{j in 1..M}:
	sum {k in 1..NU} y[k,j] = 1.0;
	
data;
param: 
	TLO 		TUP 	BLO 	BUP 	Q:= 
1	0.729961	2.11626	4.45966	397.747	250000.0
2	0.530628	1.91626	3.74950	882.353	150000.0
3	1.09024		2.47654	4.49144	833.333	180000.0
4	-0.133531	1.25276	3.14988	638.298	160000.0
5	0.0487901	1.43508	3.04452	666.667	120000.0
;

param svals:
	1	2	3	4	5	6 :=
1	7.9	2.0	5.2	4.9	6.1	4.2
2	0.7	0.8	0.9	3.4	2.1	2.5
3	0.7	2.6	1.6	3.6	3.2	2.9
4	4.7	2.3	1.6	2.7	1.2	2.5
5	1.2	3.6	2.4	4.5	1.6	2.1
;

param tvals:
        1       2       3       4       5       6 :=
1	6.4	4.7	8.3	3.9	2.1	1.2
2	6.8	6.4	6.5	4.4	2.3	3.2
3	1.0	6.3	5.4	11.9	5.7	6.2
4	3.2	3.0	3.5	3.3	2.8	3.4
5	2.1	2.5	4.2	3.6	3.7	2.2
;

solve;
display f;
display n, v, b, tl, y;
