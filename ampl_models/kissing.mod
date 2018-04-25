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

#   Source: This problem is associated to the family of Hard-Spheres 
#   problem. It belongs to the family of sphere packing problems, a 
#   class of challenging problems dating from the beginning of the 
#   17th century which is related to practical problems in Chemistry, 
#   Biology and Physics. It consists on maximizing the minimum pairwise 
#   distance between NP points on a sphere in \R^{MDIM}. 
#   This problem may be reduced to a nonconvex nonlinear optimization 
#   problem with a potentially large number of (nonoptimal) points 
#   satisfying optimality conditions. We have, thus, a class of problems 
#   indexed by the parameters MDIM and NP, that provides a suitable 
#   set of test problems for evaluating nonlinear programming codes.
#   After some algebric manipulations, we can formulate this problem as
#
#                            Minimize z
#
#                            subject to
#       
#      z \geq <x_i, x_j> for all different pair of indices i, j
#      
#                            ||x_i||^2 = 1    for all i = 1,...,NP
#
#     The goal is to find an objective value less than 0.5 (This means
#     that the NP points stored belong to the sphere and every distance
#     between two of them is greater than 1.0).
#
#     Obs: the starting point is aleatorally chosen although each 
#     variable belongs to [-1.,1.].
#
#     References:
#     [1] "Validation of an Augmented Lagrangian algorithm with a 
#          Gauss-Newton Hessian approximation using a set of 
#          Hard-Spheres problems", N. Krejic, J. M. Martinez, M. Mello 
#          and E. A. Pilotta, Tech. Report RP 29/98, IMECC-UNICAMP, 
#          Campinas, 1998.
#     [2] "Inexact-Restoration Algorithm for Constrained Optimization",
#          J. M. Martinez and E. A. Pilotta, Tech. Report, IMECC-UNICAMP, 
#          Campinas, 1998.
#     [3]  "Sphere Packings, Lattices and Groups", J. H. Conway and 
#           N. J. C. Sloane, Springer-Verlag, NY, 1988.
#
#
#     SIF input: September 29, 1998
#		 Jose Mario Martinez
#                Elvio Angel Pilotta
#

#   classification LQR2-RN-V-V

param np := 42;
param mdim := 3;
param x_init{1..np,1..mdim};

var x{i in 1..np, j in 1..mdim} := if (i <= 12) then x_init[i,j] else 0.0;

var z;

minimize f:
	z;
subject to cons1{i in 1..np-1, j in i+1..np}:
	sum {k in 1..mdim} x[i,k]*x[j,k] -z <= 0.0;
subject to cons2{i in 1..np}:
	sum {k in 1..mdim} x[i,k]^2 = 1.0;

data;
param x_init :=
1	1      -0.10890604
1	2       0.85395078
1	3      -0.45461680
2	1       0.49883922
2	2      -0.18439316
2	3      -0.04798594
3	1       0.28262888
3	2      -0.48054070
3	3       0.46715332
4	1      -0.00580106
4	2      -0.49987584
4	3      -0.44130302
5	1       0.81712540
5	2      -0.36874258
5	3      -0.68321896
6	1       0.29642426
6	2       0.82315508
6	3       0.35938150
7	1       0.09215152
7	2      -0.53564686
7	3       0.00191436
8	1       0.11700318
8	2       0.96722760
8	3      -0.14916438
9	1       0.01791524
9	2       0.17759446
9	3      -0.61875872
10	1     -0.63833630
10	2      0.80830972
10	3      0.45846734
11	1      0.28446456
11	2      0.45686938
11	3      0.16368980
12	1      0.76557382
12	2      0.16700944
12	3     -0.31647534
;
solve;
display f;
display x;
