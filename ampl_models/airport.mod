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
#   Contribution from a LANCELOT user.

#   SIF input : Rodrigo de Barros Nabholz & Maria Aparecida Diniz Ehrhardt
#               November 1994, DMA - IMECC- UNICAMP
#   Adaptation for CUTE: Ph. Toint, November 1994.

#   classification SQR2-MN-84-42

param N:=42;
param r{1..N};
param cx{1..N};
param cy{1..N};

var x{1..N} <= 10, >= -10;
var y{1..N} <= 10, >= -10;

minimize f:
	sum {i in 1..N-1, j in i+1..N} ( (x[i]-x[j])^2 + (y[i]-y[j])^2 );
subject to cons1 {i in 1..N}:
	(x[i]-cx[i])^2 + (y[i]-cy[i])^2 - r[i] <= 0;

data;
param r:=
1                  0.09 
2                  0.3
3                  0.09
4                  0.45
5                  0.5
6                  0.04
7                  0.1
8                  0.02
9                  0.02
10                 0.07
11                 0.4
12                 0.045
13                 0.05
14                 0.056
15                 0.36
16                 0.08
17                 0.07
18                 0.36
19                 0.67
20                 0.38
21                 0.37
22                 0.05
23                 0.4
24                 0.66
25                 0.05
26                 0.07
27                 0.08
28                 0.3
29                 0.31
30                 0.49
31                 0.09
32                 0.46
33                 0.12
34                 0.07
35                 0.07
36                 0.09
37                 0.05
38                 0.13
39                 0.16
40                 0.46
41                 0.25
42                 0.1;

param cx:=
1                 -6.3
2                 -7.8
3                 -9.0
4                 -7.2
5                 -5.7
6                 -1.9
7                 -3.5
8                 -0.5
9                 1.4
10                4.0
11                2.1
12                5.5
13                5.7
14                5.7
15                3.8
16                5.3
17                4.7
18                3.3
19                0.0
20                -1.0
21                -0.4
22                4.2
23                3.2
24                1.7
25                3.3
26                2.0
27                0.7
28                0.1
29                -0.1
30                -3.5
31                -4.0
32                -2.7
33                -0.5
34                -2.9
35                -1.2
36                -0.4
37                -0.1
38                -1.0
39                -1.7
40                -2.1
41                -1.8
42                0.0;

param cy:=
1                 8.0
2                 5.1
3                 2.0
4                 2.6
5                 5.5
6                 7.1
7                 5.9
8                 6.6
9                 6.1
10                5.6
11                4.9
12                4.7
13                4.3
14                3.6
15                4.1
16                3.0
17                2.4
18                3.0
19                4.7
20                3.4
21                2.3
22                1.5
23                0.5
24                -1.7
25                -2.0
26                -3.1
27                -3.5
28                -2.4
29                -1.3
30                0.0
31                -1.7
32                -2.1
33                -0.4
34                -2.9
35                -3.4
36                -4.3
37                -5.2
38                -6.5
39                -7.5
40                -6.4
41                -5.1
42                0.0;


solve;
display f;
display x, y;
