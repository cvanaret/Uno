#  NLP written by GAMS Convert at 06/20/02 11:35:39
#  
#  Equation counts
#     Total       E       G       L       N       X
#         5       1       0       4       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         4       4       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        12       6       6       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 >= 0, <= 1;
var x2 >= 0, <= 1;
var x3;

minimize obj:    x3;

subject to

e2: (0.0476666666666666 - 0.0649999999999999*x1)*exp(10*x1/(1 + 0.01*x1)) - x1
     - x3 <= 0;

e3: x1 - (0.0476666666666666 - 0.0649999999999999*x1)*exp(10*x1/(1 + 0.01*x1))
     - x3 <= 0;

e4: (0.143 + (-0.13*x1) - 0.195*x2)*exp(10*x2/(1 + 0.01*x2)) + x1 - 3*x2 - x3
     <= 0;

e5: (-(0.143 + (-0.13*x1) - 0.195*x2)*exp(10*x2/(1 + 0.01*x2))) - x1 + 3*x2
     - x3 <= 0;
