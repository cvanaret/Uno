#  NLP written by GAMS Convert at 06/20/02 11:34:28
#  
#  Equation counts
#     Total       E       G       L       N       X
#         7       5       0       2       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         7       7       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        34      24      10       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 >= -2, <= 2;
var x2 >= -2, <= 2;
var x3 >= -2, <= 2;
var x4 >= -2, <= 2;
var x5 >= -2, <= 2;
var x6;

minimize obj:    x6;

subject to

e2:    2*x1 + x2 + x3 + x4 + x5 = 6;

e3:    x1 + 2*x2 + x3 + x4 + x5 = 6;

e4:    x1 + x2 + 2*x3 + x4 + x5 = 6;

e5:    x1 + x2 + x3 + 2*x4 + x5 = 6;

e6: x1*x2*x3*x4*x5 - x6 <= 1;

e7:  - x1*x2*x3*x4*x5 - x6 <= -1;
