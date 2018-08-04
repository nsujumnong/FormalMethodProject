%% simple planar elbow manipulator

syms l1 l2 lc1 lc2 q1 q2 dq1 dq2;
syms m1 m2 g;
syms Ixx1 Ixy1 Ixz1 Iyx1 Iyy1 Iyz1 Izx1 Izy1 Izz1;
syms Ixx2 Ixy2 Ixz2 Iyx2 Iyy2 Iyz2 Izx2 Izy2 Izz2;

% Inertia matrix
I1 = [Ixx1 Ixy1 Ixz1;...
      Iyx1 Iyy1 Iyz1;...
      Izx1 Izy1 Izz1];

I2 = [Ixx2 Ixy2 Ixz2;...
      Iyx2 Iyy2 Iyz2;...
      Izx2 Izy2 Izz2];
% Jacobian matrix

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

% kinetic energy
ke = 1/2*m1*(vc1)'*vc1 + 1/2*m2*(vc2)'*vc2;



  

