#  NLP written by GAMS Convert at 06/20/02 11:34:06
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
#        14       8       6       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 := 0.25, >= 0.25, <= 1;
var x2 := 1.5, >= 1.5, <= 6.28;
var x3;

minimize obj:    x3;

subject to

e1: 0.5*sin(x1*x2) - 0.5*x1 - 0.0795774703703634*x2 - x3 <= 0;

e2: 0.920422529629637*exp(2*x1) - 5.4365636*x1 + 0.865255957591193*x2 - x3
     <= 2.5019678106022;

e3: 0.5*x1 - 0.5*sin(x1*x2) + 0.0795774703703634*x2 - x3 <= 0;

e5: 5.4365636*x1 - 0.920422529629637*exp(2*x1) - 0.865255957591193*x2 - x3
     <= -2.5019678106022;
