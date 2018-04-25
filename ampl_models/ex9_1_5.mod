#  NLP written by GAMS Convert at 06/20/02 12:12:32
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
#        33      23      10       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x2 >= 0;
var x3 >= 0;
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

minimize obj:  - x2 + 10*x3 - x4;

subject to

e2:    x2 + x3 + x5 = 1;

e3:    x2 + x4 + x6 = 1;

e4:    x3 + x4 + x7 = 1;

e5:  - x3 + x8 = 0;

e6:  - x4 + x9 = 0;

e7: x10*x5 = 0;

e8: x11*x6 = 0;

e9: x12*x7 = 0;

e10: x13*x8 = 0;

e11: x14*x9 = 0;

e12:    x10 + x12 - x13 = 1;

e13:    x11 + x12 - x14 = 1;
