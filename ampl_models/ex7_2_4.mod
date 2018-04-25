#  NLP written by GAMS Convert at 07/09/02 21:32:07
#  
#  Equation counts
#     Total       E       G       L       N       X       C
#         5       1       0       4       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         9       9       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        18       4      14       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 := 0.1, >= 0.1, <= 10;
var x2 := 0.1, >= 0.1, <= 10;
var x3 := 0.1, >= 0.1, <= 10;
var x4 := 0.1, >= 0.1, <= 10;
var x5 := 0.1, >= 0.1, <= 10;
var x6 := 0.1, >= 0.1, <= 10;
var x7 := 0.1, >= 0.1, <= 10;
var x8 := 0.1, >= 0.1, <= 10;

minimize obj: 0.4*x1^0.67/x7^0.67 + 0.4*x2^0.67/x8^0.67 - x1 - x2 + 10;

subject to

e2: 0.0588*x5*x7 + 0.1*x1 <= 1;

e3: 0.0588*x6*x8 + 0.1*x1 + 0.1*x2 <= 1;

e4: 4*x3/x5 + 2/(x3**0.71*x5) + 0.0588*x7/x3**1.3 <= 1;

e5: 4*x4/x6 + 2/(x4**0.71*x6) + 0.0588*x4**1.3*x8 <= 1;
