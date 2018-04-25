# AMPL Model by Hande Y. Benson
#
# Copyright (C) 2001 Princeton University
# All Rights Reserved
#
# Permission to use, copy, modify, and distribute this software and
# its documentation for any purpose and without fee is hereby
# granted, provided that the above copyright notice appear in all
# copies and that the copyright notice and this
# permission notice appear in all supporting documentation.                     

#   classification OUR2-AN-V-0

param n := 1000;
var x {j in 1..n} := j;

minimize obj: 
  sum {j in 1..n} ( (x[j] + x[((3*j-2) mod n)+1] + x[((7*j-3) mod n)+1])^2 +
	       4*cos(x[j] + x[((3*j-2) mod n)+1] + x[((7*j-3) mod n)+1])   );

solve;
display x;
display obj;
