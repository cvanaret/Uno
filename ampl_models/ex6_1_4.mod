#  NLP written by GAMS Convert at 06/20/02 11:52:07
#  
#  Equation counts
#     Total       E       G       L       N       X
#         5       5       0       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         7       7       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        22       4      18       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x2 := 7E-5, >= 1E-6, <= 1;
var x3 := 0.99686, >= 1E-6, <= 1;
var x4 := 0.00307, >= 1E-6, <= 1;
var x5 := 0.000474076675116379, >= 0;
var x6 := 0.997993046160137, >= 0;
var x7 := 0.0126755759820888, >= 0;

minimize obj: x2*(0.28809 + log(x2)) + x3*(log(x3) - 0.29158) + x4*(0.59336 + 
              log(x4)) + x2*(1.44805026165593*x6 + 0.989428667054834*x7) + x3*(
              1.12676386427658*x5 + 1.00363012835441*x7) + x4*(
              0.0347225450624344*x5 + 0.82681418300153*x6);

subject to

e2: x5*(x2 + 0.145002897355373*x3 + 0.989528214945409*x4) - x2 = 0;

e3: x6*(0.293701311601799*x2 + x3 + 0.646291923054068*x4) - x3 = 0;

e4: x7*(0.619143628558899*x2 + 0.239837817616513*x3 + x4) - x4 = 0;

e5:    x2 + x3 + x4 = 1;
