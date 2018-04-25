#  NLP written by GAMS Convert at 06/20/02 11:46:04
#  
#  Equation counts
#     Total       E       G       L       N       X
#         1       1       0       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         2       2       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#         2       1       1       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 := -3, >= -5, <= 5;

minimize obj: x1^6 - 15*x1^4 + 27*x1^2 + 250;


