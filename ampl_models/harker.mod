#  NLP written by GAMS Convert at 06/20/02 12:22:15
#  
#  Equation counts
#     Total       E       G       L       N       X
#         8       8       0       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#        21      21       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        61      41      20       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 >= 0;
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
var x15 := 25, >= 0;
var x16 := 25, >= 0;
var x17 := 25, >= 0;
var x18 := 25, >= 0;
var x19 := 25, >= 0;
var x20 := 25, >= 0;

minimize obj:  - (19*x15 - 0.1*x15^2 - 0.5*x18^2 - x18 - 0.005*x16^2 + 27*x16
               - 0.4*x19^2 - 2*x19 - 0.15*x17^2 + 30*x17 - 0.3*x20^2 - 1.5*x20
               - (0.166666666666667*x1^3 + x1 + 0.0666666666666667*x2^3 + 2*x2
               + 0.1*x3^3 + 3*x3 + 0.133333333333333*x4^3 + x4 + 0.1*x5^3 + 2*
              x5 + 0.0333333333333333*x6^3 + x6 + 0.0333333333333333*x7^3 + x7
               + 0.166666666666667*x8^3 + 3*x8 + 0.0666666666666667*x9^3 + 2*x9
               + 0.333333333333333*x10^3 + x10 + 0.0833333333333333*x11^3 + 2*
              x11 + 0.0666666666666667*x12^3 + 2*x12 + 0.3*x13^3 + x13 + 
              0.266666666666667*x14^3 + 3*x14));

subject to

e1:    x15 + x16 + x17 - x18 - x19 - x20 = 0;

e2:  - x1 - x2 + x5 + x8 - x15 + x18 = 0;

e3:  - x3 + x11 - x16 + x19 = 0;

e4:  - x4 + x12 - x17 + x20 = 0;

e5:    x1 - x5 - x6 - x7 + x9 + x13 = 0;

e6:    x2 + x6 - x8 - x9 - x10 + x14 = 0;

e7:    x3 + x4 + x7 + x10 - x11 - x12 - x13 - x14 = 0;
