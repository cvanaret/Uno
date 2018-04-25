#  NLP written by GAMS Convert at 06/20/02 12:11:31
#  
#  Equation counts
#     Total       E       G       L       N       X
#        13      13       0       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#        14      14       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        37      27      10       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1;
var x2;
var x4 >= 0;
var x5 >= 0;
var x6 >= 0;
var x7 >= 0;
var x8 >= 0;
var x9 >= 0;
var x10 >= 0;
var x11 >= 0;
var x12 >= 0;
var x13 >= 0;
var x14 >= 0;

minimize obj:  - 3*x1 + 2*x2 - x4;

subject to

e2:    x1 + 4*x2 - 2*x4 + x5 = 16;

e3:    3*x1 - 2*x2 + 8*x4 + x6 = 48;

e4:    x1 - 3*x2 - 2*x4 + x7 = -12;

e5:  - x1 + x8 = 0;

e6:    x1 + x9 = 4;

e7: x10*x5 = 0;

e8: x11*x6 = 0;

e9: x12*x7 = 0;

e10: x13*x8 = 0;

e11: x14*x9 = 0;

e12:    x10 + 3*x11 + x12 - x13 + x14 = 1;

e13:    2*x11 - 3*x12 = 0;
