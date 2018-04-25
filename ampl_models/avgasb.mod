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

#   classification QLR2-AN-8-10

# modified to a continuous problem

var x {1..8} >= 0.0 ,  <= 1.0 ,  := 0.5;

param a {1..8};
param b {1..7};
param c {2..8};

minimize obj:
	sum {j in 1..8} a[j]*x[j]^2 +
	sum {j in 1..7} b[j]*x[j]*x[j+1] +
	sum {j in 2..8} c[j]*x[j];

subject to con1 {j in 1..4}: x[2*j-1] + x[2*j] <= 1.0;
subject to con5 {j in 0..1}: sum {i in 1..4} x[2*i-j] <= 2.0;
subject to con7: 0 <= 2.0*x[1] + x[3] - x[7];
subject to con8: 0 <= 5.0*x[1] + 3.0*x[3] - 3.0*x[5] - x[7];
subject to con9: 0 <= x[2] - x[4] - 3.0*x[6] - 5.0*x[8];
subject to con10: 0 <= x[2] - 3.0*x[6] - 2.0*x[8];

data;
param a := 1  2.0 2  2.0 3  2.0 4  2.0 5  2.0 6  2.0 7  2.0 8  2.0 ;
param b := 1 -1.0 2 -1.0 3 -1.0 4 -1.0 5 -1.0 6 -1.0 7 -1.0 ;
param c :=        2 -2.0 3 -1.0 4 -3.0 5 -2.0 6 -4.0 7 -3.0 8 -5.0 ;

solve;

display x;
display obj;
