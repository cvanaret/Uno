#  NLP written by GAMS Convert at 06/20/02 12:14:33
#  
#  Equation counts
#     Total       E       G       L       N       X
#        13      13       0       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#        17      17       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        33      17      16       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x2 >= 0;
var x3 >= 0;
var x4 >= 0;
var x5 >= 0;
var x6 >= 0, <= 200;
var x7 >= 0, <= 200;
var x8 >= 0, <= 200;
var x9 >= 0, <= 200;
var x10 >= 0, <= 200;
var x11 >= 0, <= 200;
var x12 >= 0, <= 200;
var x13 >= 0, <= 200;
var x14 >= 0, <= 200;
var x15 >= 0, <= 200;
var x16 >= 0, <= 200;
var x17 >= 0, <= 200;

minimize obj: x2*x2 - 2*x2 + x3*x3 - 2*x3 + x4*x4 + x5*x5;

subject to

e2:  - x4 + x6 = -0.5;

e3:  - x5 + x7 = -0.5;

e4:    x4 + x8 = 1.5;

e5:    x5 + x9 = 1.5;

e6: x6*x12 = 0;

e7: x7*x13 = 0;

e8: x8*x14 = 0;

e9: x9*x15 = 0;

e10: x10*x16 = 0;

e11: x11*x17 = 0;

e12:  - 2*x2 + 2*x4 - x12 + x14 = 0;

e13:  - 2*x3 + 2*x5 - x13 + x15 = 0;
