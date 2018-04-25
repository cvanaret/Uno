#  NLP written by GAMS Convert at 06/20/02 12:13:22
#  
#  Equation counts
#     Total       E       G       L       N       X
#        12       9       0       3       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#        11      11       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        27      17      10       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x2 >= 0;
var x3 >= 0;
var x4 >= 0, <= 20;
var x5 >= 0, <= 20;
var x6 >= 0, <= 20;
var x7 >= 0, <= 20;
var x8 >= 0, <= 20;
var x9 >= 0, <= 20;
var x10 >= 0, <= 20;
var x11 >= 0, <= 20;

minimize obj: x2*x2 + (x3 - 10)*(x3 - 10);

subject to

e2:    x2 <= 15;

e3:  - x2 + x3 <= 0;

e4:  - x2 <= 0;

e5:    x2 + x3 + x4 = 20;

e6:  - x3 + x5 = 0;

e7:    x3 + x6 = 20;

e8: x4*x8 = 0;

e9: x5*x9 = 0;

e10: x6*x10 = 0;

e11: x7*x11 = 0;

e12:    2*x2 + 4*x3 + x8 - x9 + x10 = 60;
