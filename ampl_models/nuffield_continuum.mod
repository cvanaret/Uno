param c := 5;

var a >=0, <=1;
var b >=0, <=1;

maximize foo: (c+1)*((1-b)^2+(3-2*a))/4 
	        - c*((1-b)^2+(1-a)^2)/4
		- 3*(
		    a^3/24 + (1-a)*(1-(1-b)^2)/4 
#		    + (1-b)^3/3 + (1-b)^2*(b-a)/4
		    + (1-b)^2*(1-(a+b)/2)/2 - (1-b)^3/6
		    );

s.t. alessb:  a <= b;

let a:=0.03;
let b:=0.10;

display foo;

solve;

display a, b, foo;
