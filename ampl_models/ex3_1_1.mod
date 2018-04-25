#  NLP written by GAMS Convert at 06/20/02 11:42:41
#  
#  Equation counts
#     Total       E       G       L       N       X
#         7       1       0       6       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         9       9       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        21      13       8       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 := 579.19, >= 100, <= 10000;
var x2 := 1360.13, >= 1000, <= 10000;
var x3 := 5109.92, >= 1000, <= 10000;
var x4 := 182.01, >= 10, <= 1000;
var x5 := 295.6, >= 10, <= 1000;
var x6 := 217.99, >= 10, <= 1000;
var x7 := 286.4, >= 10, <= 1000;
var x8 := 395.6, >= 10, <= 1000;

minimize obj:    x1 + x2 + x3;

subject to

e2:    0.0025*x4 + 0.0025*x6 <= 1;

e3:  - 0.0025*x4 + 0.0025*x5 + 0.0025*x7 <= 1;

e4:  - 0.01*x5 + 0.01*x8 <= 1;

e5: 100*x1 - x1*x6 + 833.33252*x4 <= 83333.333;

e6: x2*x4 - x2*x7 - 1250*x4 + 1250*x5 <= 0;

e7: x3*x5 - x3*x8 - 2500*x5 <= -1250000;
