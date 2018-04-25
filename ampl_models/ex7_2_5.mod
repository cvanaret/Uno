#  NLP written by GAMS Convert at 07/05/02 23:00:42
#  
#  Equation counts
#     Total       E       G       L       N       X       C
#         7       1       0       6       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         6       6       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        30       1      29       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 := 78, >= 78, <= 102;
var x2 := 33, >= 33, <= 45;
var x3 := 27, >= 27, <= 45;
var x4 := 27, >= 27, <= 45;
var x5 := 27, >= 27, <= 45;

minimize obj: 5.3578*x3^2 + 0.8357*x1*x5 + 37.2392*x1;

subject to

e2: 2.584e-5*x3*x5 - 6.663e-5*x2*x5 - 7.34e-5*x1*x4 <= 1;

e3: 0.000853007*x2*x5 + 9.395e-5*x1*x4 - 0.00033085*x3*x5 <= 1;

e4: 1330.3294/(x2*x5) - 0.42*x1/x5 - 0.30586*x3**2/x2/x5 <= 1;

e5: 0.00024186*x2*x5 + 0.00010159*x1*x2 + 7.379e-5*x3**2 <= 1;

e6: 2275.1327/(x3*x5) - 0.2668*x1/x5 - 0.40584*x4/x5 <= 1;

e7: 0.00029955*x3*x5 + 7.992e-5*x1*x3 + 0.00012157*x3*x4 <= 1;
