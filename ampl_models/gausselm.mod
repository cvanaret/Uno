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

#   classification LOR2-AN-V-V

param n:=16;

var x{k in 1..n, i in k..n, j in k..n} := if (i==j) then 1.0 else 0.01;

minimize f:
	-x[n,n,n];

subject to conse{k in 1..n-1, i in k+1..n,j in k+1..n}:
	x[k,i,k]*x[k,k,j]/x[k,k,k] + x[k+1,i,j] - x[k,i,j] = 0;
subject to consmikk{k in 2..n-1, i in k+1..n}:
	x[k,i,k] - x[k,k,k] <= 0;
subject to consmkik{k in 2..n-1, i in k+1..n}:
	x[k,k,i] - x[k,k,k] <= 0;
subject to consmijk{k in 2..n-1, i in k+1..n, j in k+1..n}:
        x[k,i,j] - x[k,k,k] <= 0;
subject to conspikk{k in 2..n-1, i in k+1..n}:
        x[k,i,k] + x[k,k,k] >= 0;
subject to conspkik{k in 2..n-1, i in k+1..n}:
        x[k,k,i] + x[k,k,k] >= 0;
subject to conspijk{k in 2..n-1, i in k+1..n, j in k+1..n}:
        x[k,i,j] + x[k,k,k] >= 0;
subject to var_bnd{i in 1..n, j in 1..n}:
	-1.0 <= x[1,i,j] <= 1.0;
subject to var_bnd_diag{k in 1..n}:
	x[k,k,k] >= 0.0;

fix x[1,1,1] := 1.0;
solve;

