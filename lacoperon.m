%% lac operon model test
clear all
close all
clc


syms t real;
syms Ii(t) G(t);

dIi = diff(Ii,t);
dG = diff(G,t);

k2 = 4e5;
k5 = 0.6;
k_8 = 1e-5;
F1 = 0.0001;
k_2 = 0.03;
k_5 = 0.006;
k9 = 5000;
Ie = 91100;
k3 = 0.2;
k6 = 3e-6;
Ri = 0.01;
tau = 0.008;
k_3 = 60;
k7 = 3e-6;
chi = 0.002002;
mu = 2.00001;
k4 = 1;
k8 = 0.03;
eta = 0.005;

% for Ii = 1:0.2:2
%     j1 = uint8(Ii/0.2-4);
%     for G = 25:0.2:26
%         j2 = uint8(G/0.2 - 124);
%         
%         dIi(j1,j2) = -2 * k3 * Ii^2 * ( k8 * Ri * G^2 + tau ) /( k3 * Ii^2 + mu ) + 2* k_3 * F1 + ...
%             ((k5 * Ie - ( k9 + k_5 ) * Ii ) * k_2 * chi * k4 * eta *( k3 * Ii^2 + mu ) ) /...
%             ( k7 *( k2 *( k8 * Ri * G^2 + tau ) + k_2 *( k3 * Ii^2 + mu ) ) );
% 
%         dG(j1,j2) = -2* k8 * Ri * G^2 + (2* k_8 *( k8 * Ri * G^2 + tau ) ) /( k3 * Ii^2 + mu ) + ...
%             ( k9 * Ii * k_2 * chi * k4 * eta *( k3 * Ii^2 + mu ) ) / ...
%             ( k7 *( k2 *( k8 * Ri * G^2 + tau ) + k_2 *( k3 * Ii^2 + mu ) ) );
% 
%     end
% end

%%

f(1) = -2 * k3 * Ii^2 * ( k8 * Ri * G^2 + tau ) /( k3 * Ii^2 + mu ) + 2* k_3 * F1 + ...
    ((k5 * Ie - ( k9 + k_5 ) * Ii ) * k_2 * chi * k4 * eta *( k3 * Ii^2 + mu ) ) /...
    ( k7 *( k2 *( k8 * Ri * G^2 + tau ) + k_2 *( k3 * Ii^2 + mu ) ) );

f(2) = -2* k8 * Ri * G^2 + (2* k_8 *( k8 * Ri * G^2 + tau ) ) /( k3 * Ii^2 + mu ) + ...
    ( k9 * Ii * k_2 * chi * k4 * eta *( k3 * Ii^2 + mu ) ) / ...
    ( k7 *( k2 *( k8 * Ri * G^2 + tau ) + k_2 *( k3 * Ii^2 + mu ) ) );

%%
vars = [Ii, G];

odefun = odeFunction(f,vars);

initCond = [1,25];

ode15s(odefun, [0 500], initCond)