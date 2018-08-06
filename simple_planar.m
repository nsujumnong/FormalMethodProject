%% simple planar elbow manipulator
clear
clc
close

syms t real;
syms l1 l2 lc1 lc2 q1(t) q2(t) 
% syms dq1 dq2 ddq1 ddq2;
syms m1 m2 g;
syms I1 I2;
syms Ixx1 Ixy1 Ixz1 Iyx1 Iyy1 Iyz1 Izx1 Izy1 Izz1;
syms Ixx2 Ixy2 Ixz2 Iyx2 Iyy2 Iyz2 Izx2 Izy2 Izz2;

l1 = 4; l2 = 4;
lc1 = 2; lc2 = 2;
g = 9.81;
m1 = 1; m2 = 1;
I1 = 1; I2 = 1;

dq1 = diff(q1,t);
dq2 = diff(q2,t);
ddq1 = diff(q1,t,t);
ddq2 = diff(q2,t,t);

% Jacobian matrix
q = [q1; q2];
dq = [dq1; dq2];

Jvc1 = [-lc1*sin(q1) 0;...
        lc1*cos(q1) 0;...
             0      0];

vc1 = Jvc1*dq;

Jvc2 = [-l1*sin(q1)-lc2*sin(q1+q2) -lc2*sin(q1+q2);...
         l1*cos(q1)+lc2*cos(q1+q2)  lc2*cos(q1+q2);...
                   0                      0      ];
               
vc2 = Jvc2*dq;

w1 = [0 0 dq1];
w2 = [0 0 dq1+dq2];

% kinetic energy + potential energy
ke = 1/2*m1*(vc1)'*vc1 + 1/2*m2*(vc2)'*vc2;
pe1 = m1*g*lc1*sin(q1);
pe2 = m2*g*(l1*sin(q1)+lc2*sin(q1+q2));
pe = pe1 + pe2;

% Inertia matrix
d11 = m1*lc1^2 + m2*(l1^2+lc2^2+2*l1*lc2*cos(q2))+I1+I2;
d12 = m2*(lc2^2+l1*lc2*cos(q2))+I2;
d21 = d12;
d22 = m2*lc2^2+I2;

% Christoffel symbol
% c111 = 1/2*diff(d11,q1);
% c121 = 1/2*diff(d11,q2);
% c211 = c121;
% c221 = diff(d12,q2)-1/2*(diff(d22,q1));
% c112 = diff(d21,q1)-1/2*(diff(d11,q2));
% c122 = 1/2*diff(d22,q1);
% c212 = c122;
% c222 = 1/2*diff(d22,q2);

c111 = 0;
c121 = -m2*l1*lc2*sin(q2);
c211 = c121;
h = c121;
c221 = h;
c112 = -h;
c122 = 0;
c212 = c122;
c222 = 0;

% G1 = diff(pe,q1);
% G2 = diff(pe,q2);

G1 = (m1*lc1+m2*l1)*g*cos(q1)+m2*lc2*g*cos(q1+q2);
G2 = m2*lc2*g*cos(q1+q2);

tau1 = d11*ddq1 + d12*ddq2 + c121*dq1*dq2 + c211*dq2*dq1 + c221*dq2^2 + G1;
tau2 = d21*ddq1 + d22*ddq2 + c112*dq1^2 + G2;

tau = [tau1;tau2];

% Add PD controller
qd = [0;0];
dqd = [0;0];
Kp = 25;
Kp = Kp*eye(2);
Kp = Kp*(qd-q);

Kd = 15;
Kd = Kd*eye(2);
Kd = Kd*(dqd-dq);

Fv = [0.5; 0.5];
U = Kp + Kd ;

U1 = Kp + Kd - Fv.*[dq1,dq2].'; %dq positive
U2 = Kp + Kd + Fv.*[dq1,dq2].'; %dq negative
%% make ODE function
[eqs,vars1] = reduceDifferentialOrder(tau - U1 == 0,{q1 q2});
[M,F] = massMatrixForm(eqs,vars1);

f1 = M\F;

[eqs,vars2] = reduceDifferentialOrder(tau - U2 == 0,{q1 q2});
[M,F] = massMatrixForm(eqs,vars2);

f2 = M\F;


initCond = [0 0 0 0];

odefun = odeFunction(f1,vars1);
ode15s(odefun, [0 75], initCond)
title('simple dynamics 1')
xlabel('time(s)')
legend('q1','q2','dq1','dq2')

figure
odefun2 = odeFunction(f2,vars2);
ode15s(odefun, [0 75], initCond)
title('simple dynamics 2')
xlabel('time(s)')
legend('q1','q2','dq1','dq2')

