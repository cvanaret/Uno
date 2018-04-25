#  NLP written by GAMS Convert at 07/05/02 23:00:44
#  
#  Equation counts
#     Total       E       G       L       N       X       C
#         8       1       0       7       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#        11      11       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        25       5      20       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 := 0.01, >= 0.01, <= 15;
var x2 := 0.01, >= 0.01, <= 15;
var x3 := 0.01, >= 0.01, <= 15;
var x4 := 0.01, >= 0.01, <= 15;
var x5 := 0.01, >= 0.01, <= 15;
var x6 := 0.01, >= 0.01, <= 15;
var x7 := 0.01, >= 0.01, <= 15;
var x8 := 0.01, >= 0.01, <= 15;
var x9 := 0.01, >= 0.01, <= 15;
var x10 := 0.01, >= 0.01, <= 15;

minimize obj: 0.4*x4^0.67 + 0.4*x9^0.67 + x6;

subject to

e2: x3/x1/x2**1.5/x4/x5 + 5*x3*x5**1.2/x1/x2 <= 1;

e3:    0.05*x2 + 0.05*x3 <= 1;

e4: 10/x3 - x1/x3 <= 1;

e5: x8/x6/x7**1.5/x9/x10 + 5*x8*x10**1.2/x6/x7 <= 1;

e6: x7/x2 + x8/x2 <= 1;

e7: x1/x8 - x6/x8 <= 1;

e8:    x10 <= 0.1;
