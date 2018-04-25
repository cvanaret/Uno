#  NLP written by GAMS Convert at 07/05/02 23:00:43
#  
#  Equation counts
#     Total       E       G       L       N       X       C
#         3       1       0       2       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         5       5       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#         9       2       7       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 := 0.1, >= 0.1, <= 10;
var x2 := 0.1, >= 0.1, <= 10;
var x3 := 0.1, >= 0.1, <= 10;
var x4 := 0.1, >= 0.1, <= 10;

minimize obj: 0.4*x1^0.67/x3^0.67 - x1;

subject to

e2: 0.05882*x3*x4 + 0.1*x1 <= 1;

e3: 4*x2/x4 + 2/(x2**0.71*x4) + 0.05882*x3/x2**1.3 <= 1;
