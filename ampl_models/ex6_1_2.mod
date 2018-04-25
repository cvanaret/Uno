#  NLP written by GAMS Convert at 06/20/02 11:51:39
#  
#  Equation counts
#     Total       E       G       L       N       X
#         4       4       0       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#         5       5       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        13       3      10       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x2 := 0.00421, >= 1E-6, <= 1;
var x3 := 0.99579, >= 1E-6, <= 1;
var x4 := 0.0258947377097763, >= 0;
var x5 := 0.998699779997328, >= 0;

minimize obj: x2*(0.06391 + log(x2)) + x3*(log(x3) - 0.02875) + 
              0.925356626778358*x2*x5 + 0.746014540096753*x3*x4;

subject to

e2: x4*(x2 + 0.159040857374844*x3) - x2 = 0;

e3: x5*(0.307941026821595*x2 + x3) - x3 = 0;

e4:    x2 + x3 = 1;
