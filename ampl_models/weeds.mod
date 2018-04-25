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

#   Source: p. 144 in
#   J.C. Nash,
#   "Compact numerical methods for computers"
#   (2nd edition), Adam Hilger, 1990.

#   SIF input: J.C. Nash and Ph. Toint, June 1992.

#   classification SBR2-RN-3-0

param m;

set M := 1..m;

param Y{M};
param Time{M};

var B1 := 1;
var B2 := 1;
var B3 := 1, <= 3;

var expt{i in M} = exp(-1 * (B2 + ( B3 *  Time[i])));
 
minimize L2_fit:
	 sum{i in M} ((Y[i] - B1 / (1 + expt[i]) ) ^ 2);

data;
param m:=12;  
param Y := 
 1 5.308
 2 7.24
 3 9.638
 4 12.866
 5 17.069
 6 23.192
 7 31.443
 8 38.558
 9 50.156
 10 62.948
 11 75.995
 12 91.972;

param Time := 
1 1
2 2
3 3
4 4
5 5
6 6
7 7
8 8
9 9
10 10
11 11
12 12;

solve;
