#  NLP written by GAMS Convert at 06/20/02 12:15:02
#  
#  Equation counts
#     Total       E       G       L       N       X
#         6       6       0       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         7       7       0       0       0       0       0       0
#  FX     2       2       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        14       8       6       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x2 >= 0, <= 1;
var x3 >= 0;
var x4 >= 0, <= 20;
var x5 >= 0, <= 20;
var x6 >= 0, <= 0;
var x7 >= 0, <= 0;

minimize obj: 3*x3 - 4*x2*x3 + 2*x2 + 1;

subject to

e2:  - x3 + x4 = 0;

e3:    x3 + x5 = 1;

e4: x6*x4 = 0;

e5: x7*x5 = 0;

e6:    4*x2 - x6 + x7 = 1;
