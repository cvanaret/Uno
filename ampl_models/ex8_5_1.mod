#  NLP written by GAMS Convert at 07/05/02 23:00:39
#  
#  Equation counts
#     Total       E       G       L       N       X       C
#         6       5       1       0       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         7       7       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        23      11      12       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x2 := 0.6967, >= 0;
var x3 := 0.209, >= 0;
var x4 := 0.0943, >= 0;
var x5 := 0.26249, >= 0;
var x6;
var x7;

minimize obj: x2*log(x2) + x3*log(x3) + x4*log(x4) + x7/(x5 - x7) - log(x5 - x7
              ) - 2*x6/x5 + 0.430983578191493*x2 + 3.80082402249182*x3 + 
              2.92297302249182*x4;

subject to

e2: x5^3 - (1 + x7)*x5^2 + x6*x5 - x6*x7 = 0;

e3:  - (0.37943*x2*x2 + 0.75885*x2*x3 + 0.48991*x2*x4 + 0.75885*x3*x2 + 0.8836*
    x3*x3 + 0.23612*x3*x4 + 0.48991*x4*x2 + 0.23612*x4*x3 + 0.63263*x4*x4) + x6
     = 0;

e4:  - 0.14998*x2 - 0.14998*x3 - 0.14998*x4 + x7 = 0;

e5:    x2 + x3 + x4 = 1;

e6:    x5 - x7 >= 0;
