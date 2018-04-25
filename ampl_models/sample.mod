#  NLP written by GAMS Convert at 06/20/02 13:03:03
#  
#  Equation counts
#     Total       E       G       L       N       X
#         3       1       0       2       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         5       5       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        13       5       8       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 := 200, >= 100, <= 400000;
var x2 := 200, >= 100, <= 300000;
var x3 := 200, >= 100, <= 200000;
var x4 := 200, >= 100, <= 100000;

minimize obj:    x1 + x2 + x3 + x4;

subject to

e1: 4/x1 + 2.25/x2 + 1/x3 + 0.25/x4 <= 0.0401;

e2: 0.16/x1 + 0.36/x2 + 0.64/x3 + 0.64/x4 <= 0.010085;
