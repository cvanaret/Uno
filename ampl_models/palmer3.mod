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
#   M. Palmer, Edinburgh, private communication.

#   SIF input: Nick Gould, 1990.

#   classification SBR2-RN-4-0

set M := {1..23};

param X{M};
param Y{M};

var A := 1.0;
var B >= 0, := 1.0;
var C >= 0, := 1.0;
var D >= 0, := 1.0;

minimize palmer:
	sum {m in M} ( Y[m] - 
(A*(X[m]^2) + ( B / ( C + (X[m]^2)/D ) ) ) )^2;

data;  
param X:= 
1                   -1.658063
2                   -1.570796
3                   -1.396263
4                   -1.221730
5                   -1.047198
6                   -0.872665
7                   -0.766531
8                   -0.698132
9                   -0.523599
10                  -0.349066
11                  -0.174533
12                  0.0
13                  0.174533
14                  0.349066
15                  0.523599
16                  0.698132
17                  0.766531
18                  0.872665
19                  1.047198
20                  1.221730
21                  1.396263
22                  1.570796
23                  1.658063;

param Y:=
1                  64.87939
2                  50.46046
3                  28.2034
4                  13.4575
5                  4.6547
6                  0.59447
7                  0.0000
8                  0.2177
9                  2.3029
10                 5.5191
11                 8.5519
12                 9.8919
13                 8.5519
14                 5.5191
15                 2.3029
16                 0.2177
17                 0.0000
18                 0.59447
19                 4.6547
20                 13.4575
21                 28.2034
22                 50.46046
23                 64.87939;

solve;

display A,B,C,D;
