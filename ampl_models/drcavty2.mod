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
#   P.N. Brown and Y. Saad, 
#   "Hybrid Krylov Methods for Nonlinear Systems of Equations",
#   SIAM J. Sci. Stat. Comput. 11, pp. 450-481, 1990.
#   The boundary conditions have been set according to
#   I.E. Kaporin and O. Axelsson,
#   "On a class of nonlinear equation solvers based on the residual norm
#   reduction over a sequence of affine subspaces",
#   SIAM J, Sci. Comput. 16(1), 1995.

#   SIF input: Ph. Toint, Jan 1995.

#   classification NQR2-MY-V-V

param M := 100;
param H := 1/(M+2);
param RE := 1000.0;

var y{-1..M+2, -1..M+2} := 0.0;

minimize f:
	sum {i in 1..M, j in 1..M} (20*y[i,j]-8*y[i-1,j]-8*y[i+1,j]
-8*y[i,j-1]-8*y[i,j+1]+2*y[i-1,j+1]+2*y[i+1,j-1]+2*y[i-1,j-1]+2*y[i+1,j+1] +
y[i-2,j] + y[i+2,j] + y[i,j-2] + y[i,j+2] + (RE/4)*(y[i,j+1]-y[i,j-1])
*(y[i-2,j]+y[i-1,j-1]+y[i-1,j+1]-4*y[i-1,j]-4*y[i+1,j]-y[i+1,j-1] -
y[i+1,j+1] - y[i+2,j]) - (RE/4)*(y[i+1,j]-y[i-1,j])*
(y[i,j-2]+y[i-1,j-1]+y[i+1,j-1]-4*y[i,j-1]-4*y[i,j+1]-y[i-1,j+1]-y[i+1,j+1] 
- y[i,j+2]))^2;

subject to cons1{j in -1..M+2}:
	y[-1,j] = 0.0;
subject to cons2{j in -1..M+2}:
	y[0,j] = 0.0;
subject to cons3{i in 1..M}:
	y[i,-1] = 0.0;
subject to cons4{i in 1..M}:
	y[i,0] = 0.0;
subject to cons5{i in 1..M}:
	y[i,M+1] = 0.0;
subject to cons6{i in 1..M}:
	y[i,M+2] = 0.0;
subject to cons7{j in -1..M+2}:
	y[M+1,j] = -H/2;
subject to cons8{j in -1..M+2}:
	y[M+2,j] = H/2;

solve;
display f;
display y;
