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

#   Source: Hartman problem 6 in
#   L. C. W. Dixon and G. P. Szego (Eds.)
#   Towards Global Optimization
#   North Holland, 1975.
#   Paper 9, page 163.

#   SIF input: A.R. Conn May 1995

#   classification OBR2-AN-6-0

param c {1..4};
param a {1..4, 1..6};
param p {1..4, 1..6};

var x{1..6} >= 0.0, <= 1.0, := 0.2;

minimize obj: 
    - sum{i in 1..4} c[i]*exp(-sum{j in 1..6} a[i,j]*(x[j]-p[i,j])^2);

data;
param c := 1 1.0 2 1.2 3 3.0 4 3.2;
param a[tr]: 1 2 3 4 :=
 1 10.0 0.05 3.0 17.0
 2 0.05 10.0 3.5 8.0
 3 17.0 17.0 1.7 0.05
 4 3.5 0.1 10.0 10.0
 5 1.7 8.0 17.0 0.1
 6 8.0 14.0 8.0 14.0
 ;
param p[tr]: 1 2 3 4 :=
 1 0.1312 0.2329 0.2348 0.4047
 2 0.1696 0.4135 0.1451 0.8828
 3 0.5569 0.8307 0.3522 0.8732
 4 0.0124 0.3736 0.2883 0.5743
 5 0.8283 0.1004 0.3047 0.1091
 6 0.5886 0.9991 0.665 0.0381
 ;

solve;
display x;
display obj;
