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

#   Source: example 3.2.4 in
#   S. Eisenstat and H. Walker,
#   "Choosing the forcing terms in an inexact Newton method"
#   Report 6/94/75, Dept of Maths, Utah State University, 1994.

#   SIF input: Ph. Toint, July 1994.

#   classification NOR2-MN-V-V

param P := 72;
param D := -50.0;
param H := 1/(P-1);

var u{i in 1..P,j in 1..P}:=1-(i-1)*(j-1)*H^2;

minimize f: 0;
subject to cons1{i in 2..P-1, j in 2..P-2}:
	((u[i+1,j]^2+u[i-1,j]^2 + u[i,j-1]^2+u[i,j+1]^2-4*u[i,j]^2)/H^2+D*(u[i+1,j]^3-u[i-1,j]^3)/(2*H)) = 0;
subject to cons2{i in 2..P-2, j in P-1..P-1}:
	((u[i+1,j]^2+u[i-1,j]^2 + u[i,j-1]^2+u[i,j+1]^2-4*u[i,j]^2)/H^2+D*(u[i+1,j]^3- u[i-1,j]^3)/(2*H)) = 0; 
subject to cons3:
	((u[P,P-1]^2+u[P-2,P-1]^2 + u[P-1,P-2]^2+u[P-1,P]^2-4*u[P-1,P-1]^2)/H^2+D*(u[P,P-1]^3-u[P-2,P-1]^3)/(2*H)+50) = 0;

fix {j in 1..P}
	u[1,j] := 1.0;
fix {j in 1..P}
	u[P,j] := 0.0;
fix {i in 2..P-1}
	u[i,P] := 1.0;
fix {i in 2..P-1}
	u[i,1] := 0.0;

solve;
display f;
display u;
