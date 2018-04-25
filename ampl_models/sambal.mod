#  NLP written by GAMS Convert at 06/20/02 13:02:44
#  
#  Equation counts
#     Total       E       G       L       N       X
#        11      11       0       0       0       0
#  
#  Variable counts
#                 x       b       i     s1s     s2s      sc      si
#     Total    cont  binary integer    sos1    sos2   scont    sint
#        18      18       0       0       0       0       0       0
#  FX     0       0       0       0       0       0       0       0
#  
#  Nonzero counts
#     Total   const      NL     DLL
#        48      35      13       0
# 
#  Reformualtion has removed 1 variable and 1 equation


var x1 := 15;
var x2 := 3;
var x3 := 130;
var x4 := 80;
var x5;
var x6;
var x7 := 15;
var x8 := 130;
var x9 := 20;
var x10 := 25;
var x11 := 40;
var x12 := 55;
var x13 := 220;
var x14;
var x15;
var x16 := 190;
var x17 := 105;

minimize obj: 0.0666666666666667*(15 - x1)^2 + 0.333333333333333*(3 - x2)^2 + 
              0.00769230769230769*(130 - x3)^2 + 0.0125*(80 - x4)^2 + 
              0.0666666666666667*(15 - x7)^2 + 0.00769230769230769*(130 - x8)^2
               + 0.05*(20 - x9)^2 + 0.04*(25 - x10)^2 + 0.025*(40 - x11)^2 + 
              0.0181818181818182*(55 - x12)^2 + 0.00454545454545455*(220 - x13)
              ^2 + 0.00526315789473684*(190 - x16)^2 + 0.00952380952380952*(105
               - x17)^2;

subject to

e1:  - x1 - x2 - x3 - x4 + x13 = 0;

e2:  - x5 + x14 = 0;

e3:  - x6 + x15 = 0;

e4:  - x7 - x8 - x9 + x16 = 0;

e5:  - x10 - x11 - x12 + x17 = 0;

e6:  - x5 - x6 + x13 = 0;

e7:  - x1 - x7 - x10 + x14 = 0;

e8:  - x2 - x8 - x11 + x15 = 0;

e9:  - x3 - x12 + x16 = 0;

e10:  - x4 - x9 + x17 = 0;
