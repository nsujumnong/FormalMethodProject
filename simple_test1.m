
syms t real;
syms q1 q2 q3 q4;
syms q1t(t) q2t(t) q3t(t) q4t(t);

f1(1) = q3;

f1(2) = q4;

f1(3) = (6250*q1)/(3200*cos(q2)^2 - 5250) - (6250*q2 - 3875*q3 + 3875*q4 - 14715*cos(q1) + 2000*q3^2*sin(q2) + 2000*q4^2*sin(q2) + 10000*q2*cos(q2) + 6200*q4*cos(q2) + 7848*cos(q1 + q2)*cos(q2) + 4000*q3*q4*sin(q2) + 3200*q3^2*cos(q2)*sin(q2))/(3200*cos(q2)^2 - 5250);

f1(4) = (32500*q2 - 3875*q3 + 20150*q4 + 20601*cos(q1 + q2) - 14715*cos(q1) - 23544*cos(q1)*cos(q2) + 10400*q3^2*sin(q2) + 2000*q4^2*sin(q2) + 20000*q2*cos(q2) - 6200*q3*cos(q2) + 12400*q4*cos(q2) + 7848*cos(q1 + q2)*cos(q2) + 4000*q3*q4*sin(q2) + 6400*q3^2*cos(q2)*sin(q2) + 3200*q4^2*cos(q2)*sin(q2) + 6400*q3*q4*cos(q2)*sin(q2))/(3200*cos(q2)^2 - 5250) - (q1*(10000*cos(q2) + 6250))/(3200*cos(q2)^2 - 5250);

f1t = subs(f1,{q1 q2 q3 q4},{q1t q2t q3t q4t});

f2(1) = q3;

f2(2) = q4;

f2(3) = -(3625*q4 - 3625*q3 - 14715*cos(q1) - 6250*q1 + 6250*q2 + 7848*cos(q2)*cos(q1 + q2) + 5800*cos(q2)*q4 + 10000*cos(q2)*q2 + 2000*sin(q2)*q3^2 + 2000*sin(q2)*q4^2 + 3200*cos(q2)*sin(q2)*q3^2 + 4000*sin(q2)*q3*q4)/(50*(64*cos(q2)^2 - 105));

f2(4) = (20601*cos(q1 + q2) - 14715*cos(q1) - 3625*q3 + 18850*q4 - 6250*q1 + 32500*q2 + 7848*cos(q2)*cos(q1 + q2) - 5800*cos(q2)*q3 + 11600*cos(q2)*q4 - 10000*cos(q2)*q1 + 20000*cos(q2)*q2 + 10400*sin(q2)*q3^2 + 2000*sin(q2)*q4^2 - 23544*cos(q1)*cos(q2) + 6400*cos(q2)*sin(q2)*q3^2 + 3200*cos(q2)*sin(q2)*q4^2 + 4000*sin(q2)*q3*q4 + 6400*cos(q2)*sin(q2)*q3*q4)/(50*(64*cos(q2)^2 - 105));

f2t = subs(f2,{q1 q2 q3 q4},{q1t q2t q3t q4t});

%%
vars = [q1t q2t q3t q4t];
initCondition = [0 0 0 0];
odefun1 = odeFunction(f1t,vars);
[t1,y1] = ode15s(odefun1, [0 50], initCondition);
plot(t1,y1)
title('f1')

odefun2 = odeFunction(f2t,vars);
[t2,y2] = ode15s(odefun2, [0 50], initCondition);
figure
plot(t2,y2)
title('f2')