#  NLP written by GAMS Convert at 07/05/02 23:00:43
#  
#  Equation counts
#     Total       E       G       L       N       X       C
#         2       1       0       1       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         4       4       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#         6       1       5       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 := 1, >= 1, <= 100;
var x2 := 1, >= 1, <= 100;
var x3 := 1, >= 1, <= 100;

minimize obj: 0.5*x1/x2 - x1 - 5/x2;

subject to

e2: 0.01*x2/x3 + 0.0005*x1*x3 + 0.01*x1 <= 1;
