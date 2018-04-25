#
#**************************
# SET UP THE INITIAL DATA *
#**************************
#   Problem :
#   *********
#   A production planning problem in the computer industry.
#   Source:
#   L. Escudero, private communication, 1991.
#   SIF input: A.R. Conn, March 1991.
#   classification LQR2-RY-60-29
#   Constants
#   Solution
	param t := 5;
	param tm1 := -1 + (5);

	var x00101 >= 0.0;
	var x00201 >= 0.0;
	var x00301 >= 0.0;
	var x00401 >= 0.0;
	var x00102 >= 0.0;
	var x00202 >= 0.0;
	var x00302 >= 0.0;
	var x00402 >= 0.0;
	var x00103 >= 0.0;
	var x00203 >= 0.0;
	var x00303 >= 0.0;
	var x00403 >= 0.0;
	var x00104 >= 0.0;
	var x00204 >= 0.0;
	var x00304 >= 0.0;
	var x00404 >= 0.0;
	var x00105 >= 0.0;
	var x00205 >= 0.0;
	var x00305 >= 0.0;
	var x00405 >= 0.0;
	var i00101 >= 0.0;
	var y00101 >= 0.0;
	var i00201 >= 0.0;
	var y00201 >= 0.0;
	var i00301 >= 0.0;
	var y00301 >= 0.0;
	var i00401 >= 0.0;
	var y00401 >= 0.0;
	var i00102 >= 0.0;
	var y00102 >= 0.0;
	var i00202 >= 0.0;
	var y00202 >= 0.0;
	var i00302 >= 0.0;
	var y00302 >= 0.0;
	var i00402 >= 0.0;
	var y00402 >= 0.0;
	var i00103 >= 0.0;
	var y00103 >= 0.0;
	var i00203 >= 0.0;
	var y00203 >= 0.0;
	var i00303 >= 0.0;
	var y00303 >= 0.0;
	var i00403 >= 0.0;
	var y00403 >= 0.0;
	var i00104 >= 0.0;
	var y00104 >= 0.0;
	var i00204 >= 0.0;
	var y00204 >= 0.0;
	var i00304 >= 0.0;
	var y00304 >= 0.0;
	var i00404 >= 0.0;
	var y00404 >= 0.0;
	var i00105 >= 0.0;
	var y00105 >= 0.0;
	var i00205 >= 0.0;
	var y00205 >= 0.0;
	var i00305 >= 0.0;
	var y00305 >= 0.0;
	var i00405 >= 0.0;
	var y00405 >= 0.0;

minimize obj:
	i00101 + 2.0*y00101 + 2.0*i00201 + 3.0*y00201 + 3.0*i00301 + 2.0*y00301 + 
	4.0*i00401 + 5.0*y00401 + i00102 + 2.0*y00102 + 2.0*i00202 + 3.0*y00202 + 
	3.0*i00302 + 2.0*y00302 + 4.0*i00402 + 5.0*y00402 + i00103 + 2.0*y00103 + 
	2.0*i00203 + 3.0*y00203 + 3.0*i00303 + 2.0*y00303 + 4.0*i00403 + 5.0*y00403 + 
	i00104 + 2.0*y00104 + 2.0*i00204 + 3.0*y00204 + 3.0*i00304 + 2.0*y00304 + 
	4.0*i00404 + 5.0*y00404 + i00105 + 2.0*y00105 + 2.0*i00205 + 3.0*y00205 + 
	3.0*i00305 + 2.0*y00305 + 4.0*i00405 + 5.0*y00405;

subject to k01:
	0 >= x00101 + x00201 + x00301 + x00401 - 3.0;
subject to k02:
	0 >= x00102 + x00202 + x00302 + x00402 - 6.0;
subject to k03:
	0 >= x00103 + x00203 + x00303 + x00403 - 10.0;
subject to k04:
	0 >= x00104 + x00204 + x00304 + x00404 - 2000.0;
subject to k05:
	0 >= x00105 + x00205 + x00305 + x00405 - 18.0;
subject to d00101:
	x00101 - i00101 + y00101 - 1.0 = 0;
subject to d00201:
	x00201 - i00201 + y00201 - 1.0 = 0;
subject to d00301:
	x00301 - i00301 + y00301 - 1.0 = 0;
subject to d00401:
	x00401 - i00401 + y00401 - 1.0 = 0;
subject to d00102:
	x00102 + i00101 - i00102 + y00102 - 2.667 = 0;
subject to d00202:
	x00202 + i00201 - i00202 + y00202 - 1.667 = 0;
subject to d00302:
	x00302 + i00301 - i00302 + y00302 - 2.667 = 0;
subject to d00402:
	x00402 + i00401 - i00402 + y00402 - 3.333 = 0;
subject to d00103:
	x00103 + i00102 - i00103 + y00103 - 2.667 = 0;
subject to d00203:
	x00203 + i00202 - i00203 + y00203 - 2.0 = 0;
subject to d00303:
	x00303 + i00302 - i00303 + y00303 - 3.0 = 0;
subject to d00403:
	x00403 + i00402 - i00403 + y00403 - 3.0 = 0;
subject to d00104:
	x00104 + i00103 - i00104 + y00104 - 2.667 = 0;
subject to d00204:
	x00204 + i00203 - i00204 + y00204 - 2.667 = 0;
subject to d00304:
	x00304 + i00303 - i00304 + y00304 - 2.667 = 0;
subject to d00404:
	x00404 + i00403 - i00404 + y00404 - 2.667 = 0;
subject to d00105:
	x00105 + i00104 - i00105 + y00105 - 2.667 = 0;
subject to d00205:
	x00205 + i00204 - i00205 + y00205 - 2.333 = 0;
subject to d00305:
	x00305 + i00304 - i00305 + y00305 - 2.333 = 0;
subject to d00405:
	x00405 + i00404 - i00405 + y00405 - 2.333 = 0;
subject to smooth1:
	0 >= ((x00302+x00402)-x00301+x00401)^2 - 0.1 * (x00301+x00401)^2;
subject to smooth2:
	0 >= ((x00303+x00403)-x00302+x00402)^2 - 0.1 * (x00302+x00402)^2;
subject to smooth3:
	0 >= ((x00304+x00404)-x00303+x00403)^2 - 0.1 * (x00303+x00403)^2;
subject to smooth4:
	0 >= ((x00305+x00405)-x00304+x00404)^2 - 0.1 * (x00304+x00404)^2;

solve;
	display x00101;
	display x00201;
	display x00301;
	display x00401;
	display x00102;
	display x00202;
	display x00302;
	display x00402;
	display x00103;
	display x00203;
	display x00303;
	display x00403;
	display x00104;
	display x00204;
	display x00304;
	display x00404;
	display x00105;
	display x00205;
	display x00305;
	display x00405;
	display i00101;
	display y00101;
	display i00201;
	display y00201;
	display i00301;
	display y00301;
	display i00401;
	display y00401;
	display i00102;
	display y00102;
	display i00202;
	display y00202;
	display i00302;
	display y00302;
	display i00402;
	display y00402;
	display i00103;
	display y00103;
	display i00203;
	display y00203;
	display i00303;
	display y00303;
	display i00403;
	display y00403;
	display i00104;
	display y00104;
	display i00204;
	display y00204;
	display i00304;
	display y00304;
	display i00404;
	display y00404;
	display i00105;
	display y00105;
	display i00205;
	display y00205;
	display i00305;
	display y00305;
	display i00405;
	display y00405;
display obj;
