#  NLP written by GAMS Convert at 06/20/02 11:46:43
#  
#  Equation counts
#     Total       E       G       L       N       X
#         2       2       0       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         3       3       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#         5       3       2       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 := 0.7175, >= 0, <= 2;
var x2 := 1.47, >= 0, <= 3;

minimize obj: x2^2 - 7*x2 - 12*x1;

subject to

e2:  - 2*x1^4 - x2 = -2;
