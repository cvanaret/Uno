#  NLP written by GAMS Convert at 06/20/02 11:40:51
#  
#  Equation counts
#     Total       E       G       L       N       X
#         6       1       0       5       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         7       7       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        37      36       1       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 >= 0, <= 1;
var x2 := 6, >= 0;
var x3 >= 0;
var x4 := 1, >= 0, <= 1;
var x5 := 1, >= 0, <= 1;
var x6 >= 0, <= 2;

minimize obj: 6.5*x1 - 0.5*x1*x1 - x2 - 2*x3 - 3*x4 - 2*x5 - x6;

subject to

e2:    x1 + 2*x2 + 8*x3 + x4 + 3*x5 + 5*x6 <= 16;

e3:  - 8*x1 - 4*x2 - 2*x3 + 2*x4 + 4*x5 - x6 <= -1;

e4:    2*x1 + 0.5*x2 + 0.2*x3 - 3*x4 - x5 - 4*x6 <= 24;

e5:    0.2*x1 + 2*x2 + 0.1*x3 - 4*x4 + 2*x5 + 2*x6 <= 12;

e6:  - 0.1*x1 - 0.5*x2 + 2*x3 + 5*x4 - 5*x5 + 3*x6 <= 3;
