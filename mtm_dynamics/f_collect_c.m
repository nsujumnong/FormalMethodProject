%   /*
%   V1.2
%   */

syms t real;
syms q2t(t);
syms Dq1tt1(t) Dq2tt1(t) Dq3tt1(t);
  t3 = q2t(t);
  t7 = atan(4.575003934450549);
  t8 = t3+t7;
  t2 = cos(t8);
  t4 = Dq1tt1(t);
  t5 = Dq2tt1(t);
  t6 = t5*t5;
  t9 = cos(t3);
  t10 = 1.102114685206562E16;
  t11 = t3*2.0;
  t12 = t4*t4;
  t13 = sin(t3);
  t14 = t2*t2;
  t15 = t14*1.897901217731186E36;
  t16 = 7.601698647277576E17;
  t17 = atan(7.986526751911008E-1);
  t18 = t11+t17;
  t19 = cos(t18);
  t20 = t16*t19*2.305948830758444E20;
  t21 = t15+t20-2.836460967663414E38;
  t22 = 1.0/t21;
  t23 = 8.453429026518819E15;
  t24 = atan(1.931646285855638);
  t25 = -t3+t24;
  t26 = cos(t25);
  t27 = cos(t11);
  t28 = sin(t11);
  t29 = t9*t9;
  
  A0(1) = t4;
  A0(2) = t5;
  A0(3) = t22*(t4*2.442755675041219E36+t4*t5*1.093910323337682E38-t6*t9*9.698444859689275E36-t6*t13*2.119876834784415E36+t2*t5*t10*1.247357707171888E20-t4*t5*t29*2.187820646675363E38-t4*t5*t9*t13*2.739389367414144E38+t2*t10*t12*t27*4.74386208725146E20+t2*t10*t12*t28*5.939831211504242E20-t2*t10*t23*t26*4.905E6)*-2.0;
  A0(4) = t22*(t5*-7.855758453005507E37-t12*t27*2.987646164171462E38-t12*t28*3.740857893522464E38+t23*t26*3.089129524790131E24-t2*t4*t10*8.474622307166383E19+t5*t16*t19*6.386471460720066E19-t2*t4*t5*t10*3.795089669801168E21+t2*t6*t9*t10*3.364669581674762E20+t2*t6*t10*t13*7.354462706224653E19+t12*t16*t19*t27*2.428857388672748E20+t12*t16*t19*t28*3.041193580290172E20-t16*t19*t23*t26*2.51136E6+t2*t4*t5*t10*t29*7.590179339602336E21+t2*t4*t5*t9*t10*t13*9.503729938406787E21)*(1.0/2.0);
