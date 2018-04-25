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

#   classification QLR2-AN-V-V

param nx := 100;
param ny := 100;

var x{1..nx,0..ny+1} >= 0;
var y{0..nx+1,1..ny} >= 0;

minimize f:
	(sum {i in 1..nx-1, j in 1..ny-1}  (x[i,j] - 1)^2 #ox(i,j)
	+ sum {i in 1..nx-1, j in 1..ny-1} (y[i,j] - 1)^2 #oy(i,j)
	+ sum {i in 1..nx-1} (x[i,ny] - 1)^2 #ox(i,ny)
	+ sum {j in 1..ny-1} (y[nx,j] - 1)^2 #oy(nx,j)
	+ sum {i in 1..nx} (x[i,0] - 1)^2 #ox(i,0)
	+ sum {i in 1..nx} (x[i,ny+1] - 1)^2 #ox(i,ny+1)
	+ sum {j in 1..ny} (y[0,j] - 1)^2 #oy(0,j)
	+ sum {j in 1..ny} (y[nx+1,j] - 1)^2 #oy(nx+1,j)
	)/2;

subject to v1{i in 2..nx-1, j in 2..ny-1}:
	(x[i,j] - x[i-1,j]) + (y[i,j] - y[i,j-1]) - 1 = 0;
subject to v2{i in 2..nx-1}:
	x[i,0] + (x[i,1] - x[i-1,1]) + y[i,1] - 1 = 0;
subject to v3{i in 2..nx-1}:
        x[i,ny+1] + (x[i,ny] - x[i-1,ny]) - y[i,ny-1] - 1 = 0;
subject to v4{j in 2..ny-1}:    
	y[0,j] + (y[1,j] - y[1,j-1]) + x[1,j] - 1 = 0; 
subject to v5{j in 2..ny-1}:            
        y[nx+1,j] + (y[nx,j] - y[nx,j-1]) - x[nx-1,j] - 1 = 0;
subject to bound1{i in 1..nx-1}:
	x[i,ny] >= 1.0;
subject to bound2{j in 1..ny-1}:
	y[nx,j] >= 1.0;

solve;
