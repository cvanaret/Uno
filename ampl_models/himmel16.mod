#  NLP written by GAMS Convert at 06/20/02 12:23:09
#  
#  Equation counts
#     Total       E       G       L       N       X
#        22       7       0      15       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#        19      19       0       0       0       0       0       0
#  FX     3       3       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        97      13      84       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 >= 0, <= 0;
var x2 := 0.5;
var x3 := 0.5;
var x4 := 0.5;
var x5;
var x6;
var x7 >= 0, <= 0;
var x8 >= 0, <= 0;
var x9 := 0.4;
var x10 := 0.8;
var x11 := 0.8;
var x12 := 0.4;
var x13;
var x14;
var x15;
var x16;
var x17;
var x18;

minimize obj:  - x13 - x14 - x15 - x16 - x17 - x18;

subject to

e1: (x1 - x2)^2 + (x7 - x8)^2 <= 1;

e2: (x1 - x3)^2 + (x7 - x9)^2 <= 1;

e3: (x1 - x4)^2 + (x7 - x10)^2 <= 1;

e4: (x1 - x5)^2 + (x7 - x11)^2 <= 1;

e5: (x1 - x6)^2 + (x7 - x12)^2 <= 1;

e6: (x2 - x3)^2 + (x8 - x9)^2 <= 1;

e7: (x2 - x4)^2 + (x8 - x10)^2 <= 1;

e8: (x2 - x5)^2 + (x8 - x11)^2 <= 1;

e9: (x2 - x6)^2 + (x8 - x12)^2 <= 1;

e10: (x3 - x4)^2 + (x9 - x10)^2 <= 1;

e11: (x3 - x5)^2 + (x9 - x11)^2 <= 1;

e12: (x3 - x6)^2 + (x9 - x12)^2 <= 1;

e13: (x4 - x5)^2 + (x10 - x11)^2 <= 1;

e14: (x4 - x6)^2 + (x10 - x12)^2 <= 1;

e15: (x5 - x6)^2 + (x11 - x12)^2 <= 1;

e17:  - 0.5*(x1*x8 - x7*x2) + x13 = 0;

e18:  - 0.5*(x2*x9 - x8*x3) + x14 = 0;

e19:  - 0.5*(x3*x10 - x9*x4) + x15 = 0;

e20:  - 0.5*(x4*x11 - x10*x5) + x16 = 0;

e21:  - 0.5*(x5*x12 - x11*x6) + x17 = 0;

e22:  - 0.5*(x6*x7 - x12*x1) + x18 = 0;
