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

#   classification OLR2-RN-49-15

var SR >=10,<=150,:=2.7452D+01;
var LR >=0,<=10,:=1.5000D+00;
var PK >=0,<=10,:=1.0000D+01;
var EF >=0,<=5,:=0.0000D+00;
var SX >=7,<=120,:=1.9217D+01;
var LX >=1.5,<=8,:=1.5000D+00;
var SD >=2,<=20,:=3.5688D+00;
var SK >=2,<=30,:=4.0696D+00;
var ST >=30,<=500,:=3.4315D+01;
var SF >=20,<=200,:=8.8025D+01;
var LF >=0.01,<=20,:=5.1306D+00;
var AM >=0,<=10,:=0.0000D+00;
var CA >=-0.2,<=-0.001,:=-1.4809D-01;
var CB >=0.1,<=2,:=7.5980D-01;
var SO >=0,<=1,:=0.0000D+00;
var SS >=0,<=2,:=0.0000D+00;
var IMPDER >= 100,<=1000,:=1.1470D+02;
var IMPK >=500,<=5000,:=5.0000D+02;
var IMPFUS >=500,<=5000,:=1.7605D+03;
var QI >=1000,<=20000,:=2.3256D+03;
var PT >=2,<=30,:=5.6788D+00;
var MV >=2000,<=20000,:=1.4197D+04;
var MC >=3000,<=30000,:=1.2589D+04;
var MD >=5000,<=50000,:=2.8394D+04;
var PD >=0.2,<=0.8,:=2.0000D-01;
var NS >=1,<=5,:=1.0000D+00;
var VS >=0,<=20,:=0.0000D+00;
var CR >=100,<=400,:=1.0000D+02;
var PM >=4,<=15,:=1.5000D+01;
var DV >=0,<=10,:=0.0000D+00;
var MZ >=500,<=10000,:=5.0000D+02;
var VN >=10,<=50,:=1.0000D+01;
var QV >=250,<=5000,:=8.1490D+02;
var QF >=750,<=15000,:=3.1405D+03;
var IMPTRAIN >= 250,<=3000,:=1.9450D+03;
var IMPMOT >=10,<=5000,:=1.9085D+02;
var IMPNMOT >=35,<=70,:=3.5000D+01;
var IMPPET >=100,<=3000,:=1.0000D+02;
var IMPPIL >=200,<=400,:=2.0000D+02;
var IMPCAN >=120,<=240,:=1.2000D+02;
var IMPSNA >=700,<=1900,:=7.0000D+02;
var MS >=100,<=1000,:=1.0000D+03;
var EL >=2,<=20,:=4.9367D+00;
var DE >=0,<=1,:=0.0000D+00;
var DS >=0,<=2,:=0.0000D+00;
var IMPVOIL >= 500,<=5000,:=5.0000D+03;
var NM >=1,<=2,:=1;
var NP >=1,<=2,:=1;
var NG >=1,<=2,:=1;

minimize f:
	(SK - 0.01*PK*SR)^2
	+ (CA - (SS-SO-CB*LF)/(LF^2) )^2
	+ (-2*AM+SO+SS + 0.01*EF/LF)^2
	+ (AM - 0.025*SO*CB^2/CA)^2
	+ (IMPDER - 27.5*SD - 1.3*SD^2)^2
	+ (IMPK - 70*SK + 8.6*SK^2)^2
	+ (QI - 1000 + MV^2/24000)^2
	+ (1000*PT - MD*PD)^2
	+ (VN + VS +QF/790 + 2 - MZ/CR +DV*PT)^2
	+ (IMPMOT - 1000*PT/(PM+20) - 12*sqrt(PT))^2
	+ (ST - 1.25*SR*NM)^2
	+ (SR - MD/MS)^2
	+ (QV - 2.4*SX*sqrt(SX)*EL/sqrt(LX))^2
	+ (SO - 0.785*DE^2*PT)^2
	+ (SS - 0.785*DS^2*PT)^2
	+ (CB - 2*(VN-CA*LF^3)/(LF^2*(3-SO*LF)))^2
	+ (IMPVOIL - 1.15*SX*(15+0.15*SX)*(8+(MC*LX/(50*SR*EL))^1.5))^2
;
subject to cons1:
	SD-0.13*SR = 0;
subject to cons2:
	SX-0.7*SR = 0;
subject to cons3:
	LX-LR = 0;
subject to cons5:
	SF - ST - 2*SD - 2*SX - 2*SK = 0;
subject to cons11:
	IMPFUS - 20*SF = 0;
subject to cons12:
	MD - 2*MV = 0;
subject to cons15:
	QF - QI - QV = 0;
subject to cons17:
	IMPTRAIN - 0.137*MV = 0;
subject to cons19:
	IMPNMOT - 35*NM = 0;
subject to cons20:
	IMPPET - 0.043*QI = 0;
subject to cons21:
	IMPPIL - 200*NP = 0;
subject to cons22:
	IMPCAN - 120*NG = 0;
subject to cons23:
	IMPSNA - 300*NS -400 = 0;
subject to cons24:
	MC - MV + 95*NP + 70*NG + 660*NM + 0.5*QI -380= 0;
subject to cons25:
	MZ - IMPTRAIN + IMPNMOT + IMPPET + IMPPIL + IMPCAN + IMPSNA
	+ 290=0;

solve;
display f;
