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

#   Source: problem from Kortanek and No
#   The problem is really a semi-infinite QP
#   to appear in SIAM J. Optimization.

#   The structure is the following :


#     min  "Sum"{ Xj^2/(2j) + Xj/j  ;  j=1,...,n }   subject to

#      "Sum"{ t^(j-1)#Xj } ; j=1,...,n  >=  b(t) for all t in [0 1].


# Four examples are considered for n = 20, corresponding to the RHS
# function, b(t) : sin(t), 1/(2-t), exp(t), and tan(t).

# The interval [0 1] is dicretized via steps of 1/1000

#   SIF input: A.R. Conn, May 1993

#   classification QLR2-AN-20-1001

param n := 20;
param m := 1000;

var x {1..n} := 2;

minimize obj: sum {j in 1..n} ( x[j]^2/(2*j) + x[j]/j);

subject to c{i in 0..m}: sum {j in 1..n} (i/m)^(j-1)*x[j] >= sin(i/m);

solve;
display x;
display obj;
