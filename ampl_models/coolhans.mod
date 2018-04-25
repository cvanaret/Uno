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
#   S. Ceria, private communication, 1995.

#   SIF input: Ph. Toint, Feb 1995.

#   classification NQR2-RN-9-9

set N;

param A{N,N};
param B{N,N};
param C{N,N};

var X{N,N} := 0.0;

var AXX{i in N, j in N} = sum{m in N} ( (sum{k in N} (A[i,k]*X[k,m])) *X[m,j] );
var BX{i in N, j in N} = sum{k in N} (B[i,k]*X[k,j]);

minimize f: 0;
subject to matrix{i in N, j in N}:
	(AXX[i,j] + BX[i,j] + C[i,j]) = 0;

data;
set N := 1 2 3;

param A :=
1 1 0
2 1 0.00000013725
3 1 0
1 2 0
2 2 937.62
3 2 0
1 3 0
2 3 -42.207
3 3 0;

param B :=
1 1 0.0060893
2 1 0.00000013880
3 1 -0.00000013877
1 2 -44.292
2 2 -1886.0
3 2 42.362
1 3 2.0011
2 3 42.362
3 3 -2.0705;

param C :=
1 1 0
2 1 0
3 1 0
1 2 44.792
2 2 948.21
3 2 -42.684
1 3 0
2 3 0
3 3 0;


solve;

display X;
