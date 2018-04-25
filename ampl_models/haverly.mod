#  NLP written by GAMS Convert at 06/20/02 12:22:29
#  
#  Equation counts
#     Total       E       G       L       N       X
#        10       8       0       2       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#        13      13       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        34      27       7       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 >= 0;
var x2 >= 0;
var x3 >= 0;
var x4 >= 0;
var x5 >= 0;
var x6 >= 0, <= 100;
var x7 >= 0, <= 200;
var x8 := 1, >= 0;
var x9 := 1, >= 0;
var x10 := 1, >= 0;
var x11 := 1, >= 0;
var x12 := 1, >= 0;

minimize obj:    x1 - x2;

subject to

e1:    x1 - 6*x3 - 16*x4 - 10*x5 = 0;

e2:    x2 - 9*x6 - 15*x7 = 0;

e3:    x6 - x8 - x10 = 0;

e4:    x7 - x9 - x11 = 0;

e5:    x3 + x4 - x10 - x11 = 0;

e6:    x5 - x8 - x9 = 0;

e7: x12*(x10 + x11) - 3*x3 - x4 = 0;

e8: x12*x10 - 2.5*x10 - 0.5*x8 <= 0;

e9: x12*x11 - 1.5*x11 + 0.5*x9 <= 0;
