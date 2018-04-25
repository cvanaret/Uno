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

#   classification OXR2-MY-64-0

param p:=7;
var x{1..p+1,1..p+1};
var A{i in 1..p, j in 1..p} = (x[i,j]-x[i+1,j+1])^2;
var B{i in 1..p, j in 1..p} = (x[i,j+1]-x[i+1,j])^2;
minimize f:
	sum {i in 1..p} sum {j in 1..p} sqrt( (1+0.5*A[i,j]*p^2+0.5*B[i,j]*p^2) )/p^2;
subject to cons1{j in 1..p+1}:
	x[1,j] = 1.0;	
subject to cons2{j in 1..p+1}:
	x[p+1,j] = 1.0;	
subject to cons3{j in 1..p+1}:
	x[j,1] = 1.0;	
subject to cons4{j in 1..p+1}:
	x[j,p+1] = 1.0;	

solve; display f; display x;
