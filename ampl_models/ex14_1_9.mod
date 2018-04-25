#  NLP written by GAMS Convert at 06/20/02 11:36:11
#  
#  Equation counts
#     Total       E       G       L       N       X
#         3       1       0       2       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         3       3       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#         6       4       2       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 := 100, >= 100, <= 1000;
var x2;

minimize obj:    x2;

subject to

e2: 4510067.11409396*x1*exp(-7548.11926028431/x1) + 0.00335570469798658*x1 - 
    2020510067.11409*exp(-7548.11926028431/x1) - x2 <= 1;

e3: (-4510067.11409396*x1*exp(-7548.11926028431/x1)) - 0.00335570469798658*x1
     + 2020510067.11409*exp(-7548.11926028431/x1) - x2 <= -1;
