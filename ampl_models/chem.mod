#  NLP written by GAMS Convert at 06/20/02 11:28:14
#  
#  Equation counts
#     Total       E       G       L       N       X
#         5       5       0       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#        12      12       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        37      26      11       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 := 0.001, >= 0.001;
var x2 := 0.001, >= 0.001;
var x3 := 0.001, >= 0.001;
var x4 := 0.001, >= 0.001;
var x5 := 0.001, >= 0.001;
var x6 := 0.001, >= 0.001;
var x7 := 0.001, >= 0.001;
var x8 := 0.001, >= 0.001;
var x9 := 0.001, >= 0.001;
var x10 := 0.001, >= 0.001;
var x11 := 0.01, >= 0.01;

minimize obj: x1*(log(x1/x11) - 6.05576803624071) + x2*(log(x2/x11) - 
              17.1307680362407) + x3*(log(x3/x11) - 34.0207680362407) + x4*(
              log(x4/x11) - 5.88076803624071) + x5*(log(x5/x11) - 
              24.6877680362407) + x6*(log(x6/x11) - 14.9527680362407) + x7*(
              log(x7/x11) - 24.0667680362407) + x8*(log(x8/x11) - 
              10.6747680362407) + x9*(log(x9/x11) - 26.6287680362407) + x10*(
              log(x10/x11) - 22.1447680362407);

subject to

e1:    x1 + 2*x2 + 2*x3 + x6 + x10 = 2;

e2:    x4 + 2*x5 + x6 + x7 = 1;

e3:    x3 + x7 + x8 + 2*x9 + x10 = 1;

e5:  - x1 - x2 - x3 - x4 - x5 - x6 - x7 - x8 - x9 - x10 + x11 = 0;
