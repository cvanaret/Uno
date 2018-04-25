#  NLP written by GAMS Convert at 06/20/02 11:57:25
#  
#  Equation counts
#     Total       E       G       L       N       X
#         1       1       0       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         3       3       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#         3       1       2       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 >= -1, <= 2;
var x2 >= -1, <= 1;

minimize obj: cos(x1)*sin(x2) - x1/(1 + x2^2);


