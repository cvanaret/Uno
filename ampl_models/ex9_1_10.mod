#  NLP written by GAMS Convert at 06/20/02 12:11:51
#  
#  Equation counts
#     Total       E       G       L       N       X
#        13      12       0       1       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#        15      15       0       0       0       0       0       0
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
var x15 >= 0;

minimize obj:  - 2*x2 + x3 + 0.5*x4;

subject to

e2:    x2 + x3 <= 2;

e3:  - 2*x2 + x4 - x5 + x6 = -2.5;

e4:    x2 - 3*x3 + x5 + x7 = 2;

e5:  - x4 + x8 = 0;

e6:  - x5 + x9 = 0;

e7: x11*x6 = 0;

e8: x12*x7 = 0;

e9: x13*x8 = 0;

e10: x14*x9 = 0;

e11: x15*x10 = 0;

e12:    x11 - x13 = 4;

e13:    x11 + x12 - x14 = -1;
