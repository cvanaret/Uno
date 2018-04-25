#  NLP written by GAMS Convert at 07/05/02 23:04:19
#  
#  Equation counts
#     Total       E       G       L       N       X       C
#        10       1       0       9       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#        12      12       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        30       4      26       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 := 0.01, >= 0.01, <= 10;
var x2 := 0.01, >= 0.01, <= 10;
var x3 := 0.01, >= 0.01, <= 10;
var x4 := 0.01, >= 0.01, <= 10;
var x5 := 0.01, >= 0.01, <= 10;
var x6 := 0.01, >= 0.01, <= 10;
var x7 := 0.01, >= 0.01, <= 10;
var x8 := 0.01, >= 0.01, <= 10;
var x9 := 0.01, >= 0.01, <= 10;
var x10 := 0.01, >= 0.01, <= 10;
var x11 := 0.01, >= 0.01, <= 10;

minimize obj: 1/x3;

subject to

e2: x7*x10 + 0.1*x10 <= 1;

e3: 10*x1*x4 + 10*x1*x4*x7**2 <= 1;

e4: 1/x4 - 100*x7*x10 <= 1;

e5: x10/x11 - 10*x8 <= 1;

e6: x2*x5/x1 - x2*x5*x8**2/x1 <= 1;

e7: 1/x5 - 10*x8*x11/x1 <= 1;

e8:  - 10*x9 + 10*x11 <= 1;

e9: x3*x6/x2 - x3*x6*x9**2/x2 <= 1;

e10: 1/x6 - x9/x2 <= 1;
