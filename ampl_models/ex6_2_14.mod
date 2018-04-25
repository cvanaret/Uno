#  NLP written by GAMS Convert at 06/20/02 11:53:36
#  
#  Equation counts
#     Total       E       G       L       N       X
#         3       3       0       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         5       5       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#         9       5       4       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x2 := 0.0583, >= 1E-7, <= 0.5;
var x3 := 0.4417, >= 1E-7, <= 0.5;
var x4 := 0.408, >= 1E-7, <= 0.5;
var x5 := 0.092, >= 1E-7, <= 0.5;

minimize obj: x2*(log(x2/(x2 + x4)) + log(x2/(x2 + 0.095173*x4))) + x4*(log(x4/
              (x2 + x4)) + log(x4/(0.30384*x2 + x4))) + (x2 + 2.6738*x4)*log(x2
               + 2.6738*x4) + (0.374*x2 + x4)*log(0.374*x2 + x4) + 2.6738*x4*
              log(x4/(x2 + 2.6738*x4)) + 0.374*x2*log(x2/(0.374*x2 + x4)) + x3*
              (log(x3/(x3 + x5)) + log(x3/(x3 + 0.095173*x5))) + x5*(log(x5/(x3
               + x5)) + log(x5/(0.30384*x3 + x5))) + (x3 + 2.6738*x5)*log(x3 + 
              2.6738*x5) + (0.374*x3 + x5)*log(0.374*x3 + x5) + 2.6738*x5*log(
              x5/(x3 + 2.6738*x5)) + 0.374*x3*log(x3/(0.374*x3 + x5)) - 3.6838*
              x2*log(x2) - 1.59549*x4*log(x4) - 3.6838*x3*log(x3) - 1.59549*x5*
              log(x5);

subject to

e2:    x2 + x3 = 0.5;

e3:    x4 + x5 = 0.5;
