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

#   classification OOR2-MY-6-100

param np := 50;
param x{1..np};
param y{1..np};

var v1 := -40.0;
var w1 := 5.0;
var d >= 1e-8, := 1.0;
var a >= 1.0, := 2.0;
var t >= 0.0, <= 6.2831852, := 1.5;
var r >= 0.39, := 0.75;

minimize f:
	(d+r)^2*acos(-( (a*d)^2 - (a*d+r)^2 +(d+r)^2)/(2*(d+r)*a*d))
	-(a*d+r)^2*acos(( (a*d)^2+ (a*d+r)^2 -(d+r)^2)/(2*(a*d+r)*a*d))
	+(d+r)*a*d*sin(acos(-( (a*d)^2 - (a*d+r)^2 +(d+r)^2)/(2*(d+r)*a*d)));
subject to cons1{i in 1..np}:
	(v1+a*d*cos(t)-x[i])^2 + (w1+a*d*sin(t)-y[i])^2 - (d+r)^2<= 0.0;
subject to cons2{i in 1..np}:
	(v1-x[i])^2 + (w1-y[i])^2 - (a*d+r)^2 >= 0.0;
data;
param: x y :=
1                  0.514
0.176

2                  0.948
0.172

3                  0.702
0.226

4                  0.495
0.125

5                  0.823
0.152

6                  0.625
0.315

7                  0.347
0.917

8                  0.520
0.401

9                  0.607
0.785

10                 0.758
0.582

11                 0.200
0.827

12                 0.416
0.464

13                 0.979
0.126

14                 0.213
0.958

15                 0.737
0.409

16                 0.957
0.028

17                 0.319
0.757

18                 0.572
0.119

19                 0.570
0.252

20                 0.496
0.237

21                 0.477
0.406

22                 0.873
0.427

23                 0.522
0.697

24                 0.773
0.245

25                 0.887
0.037

26                 0.651
0.399

27                 0.676
0.733

28                 0.938
0.233

29                 0.779
0.431

30                 0.750
0.208

31                 0.803
0.219

32                 0.563
0.716

33                 0.653
0.604

34                 0.790
0.079

35                 0.246
0.945

36                 0.477
0.800

37                 0.744
0.381

38                 0.480
0.527

39                 0.446
0.705

40                 0.095
0.963

41                 0.551
0.740

42                 0.579
0.638

43                 0.782
0.188

44                 0.684
0.293

45                 0.565
0.418

46                 0.566
0.488

47                 0.607
0.416

48                 0.036
0.977

49                 0.647
0.350

50                 0.553
0.358
;

solve;
display f;
display v1, w1, a, d, r, t;
