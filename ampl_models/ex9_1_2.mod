#  NLP written by GAMS Convert at 06/20/02 12:12:02
#  
#  Equation counts
#     Total       E       G       L       N       X
#        10      10       0       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#        11      11       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        26      18       8       0
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

minimize obj:  - x2 - 3*x3;

subject to

e2:  - x2 + x3 + x4 = 3;

e3:    x2 + 2*x3 + x5 = 12;

e4:    4*x2 - x3 + x6 = 12;

e5:  - x3 + x7 = 0;

e6: x8*x4 = 0;

e7: x9*x5 = 0;

e8: x10*x6 = 0;

e9: x11*x7 = 0;

e10:    x8 + 2*x9 - x10 - x11 = -1;
