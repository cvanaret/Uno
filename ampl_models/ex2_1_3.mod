#  NLP written by GAMS Convert at 06/20/02 11:40:30
#  
#  Equation counts
#     Total       E       G       L       N       X
#        10       1       0       9       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#        14      14       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        41      37       4       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 := 1, >= 0, <= 1;
var x2 := 1, >= 0, <= 1;
var x3 := 1, >= 0, <= 1;
var x4 := 1, >= 0, <= 1;
var x5 := 1, >= 0, <= 1;
var x6 := 1, >= 0, <= 1;
var x7 := 1, >= 0, <= 1;
var x8 := 1, >= 0, <= 1;
var x9 := 1, >= 0, <= 1;
var x10 := 3, >= 0;
var x11 := 3, >= 0;
var x12 := 3, >= 0;
var x13 := 1, >= 0, <= 1;

minimize obj: 5*x1 - 0.5*(10*x1*x1 + 10*x2*x2 + 10*x3*x3 + 10*x4*x4) + 5*x2 + 5
              *x3 + 5*x4 - x5 - x6 - x7 - x8 - x9 - x10 - x11 - x12 - x13;

subject to

e2:    2*x1 + 2*x2 + x10 + x11 <= 10;

e3:    2*x1 + 2*x3 + x10 + x12 <= 10;

e4:    2*x2 + 2*x3 + x11 + x12 <= 10;

e5:  - 8*x1 + x10 <= 0;

e6:  - 8*x2 + x11 <= 0;

e7:  - 8*x3 + x12 <= 0;

e8:  - 2*x4 - x5 + x10 <= 0;

e9:  - 2*x6 - x7 + x11 <= 0;

e10:  - 2*x8 - x9 + x12 <= 0;
