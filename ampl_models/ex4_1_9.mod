#  NLP written by GAMS Convert at 06/20/02 11:46:55
#  
#  Equation counts
#     Total       E       G       L       N       X
#         3       1       0       2       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         3       3       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#         7       5       2       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 := 2.3295, >= 0, <= 3;
var x2 := 3.17846, >= 0, <= 4;

minimize obj:  - x1 - x2;

subject to

e2: 8*x1^3 - 2*x1^4 - 8*x1^2 + x2 <= 2;

e3: 32*x1^3 - 4*x1^4 - 88*x1^2 + 96*x1 + x2 <= 36;
