% This is a modification of Yan's code
% combine with a snippet of ODE function generating code from Radian's
% impedance control code

syms q1 q2 q3 q4 q5 q6 q7 real;
syms dq1 dq2 dq3 dq4 dq5 dq6 dq7 real;
syms ddq1 ddq2 ddq3 ddq4 ddq5 ddq6 ddq7 real;

%%
syms t real;
syms q1t(t) q2t(t) q3t(t) q4t(t) q5t(t) q6t(t) q7t(t);
syms dq1t(t) dq2t(t) dq3t(t) dq4t(t) dq5t(t) dq6t(t) dq7t(t);
syms ddq1t(t) ddq2t(t) ddq3t(t) ddq4t(t) ddq5t(t) ddq6t(t) ddq7t(t);
%%

% dq1t = diff(q1t,t);
% dq2t = diff(q2t,t);
% dq3t = diff(q3t,t);
% dq4t = diff(q4t,t);
% dq5t = diff(q5t,t);
% dq6t = diff(q6t,t);
% dq7t = diff(q7t,t);
% 
% ddq1t = diff(dq1t,t);
% ddq2t = diff(dq2t,t);
% ddq3t = diff(dq3t,t);
% ddq4t = diff(dq4t,t);
% ddq5t = diff(dq5t,t);
% ddq6t = diff(dq6t,t);
% ddq7t = diff(dq7t,t);


%%
% original one
A = h_b;

% vpa(A,4)
% A = [
% [ ddq2*sin(q2)*sin(q3) - 1.0*ddq3*cos(q2)*cos(q3) - 1.0*ddq2*cos(q2)*cos(q3) + ddq3*sin(q2)*sin(q3) + dq2^2*cos(q2)*sin(q3) + dq2^2*cos(q3)*sin(q2) + dq3^2*cos(q2)*sin(q3) + dq3^2*cos(q3)*sin(q2) + 2.0*dq2*dq3*cos(q2)*sin(q3) + 2.0*dq2*dq3*cos(q3)*sin(q2), 2.0*dq1*dq2 + 2.0*dq1*dq3 - 4.0*dq1*dq2*cos(q2)^2 - 4.0*dq1*dq2*cos(q3)^2 - 4.0*dq1*dq3*cos(q2)^2 - 4.0*dq1*dq3*cos(q3)^2 - 2.0*ddq1*cos(q2)*sin(q2) - 2.0*ddq1*cos(q3)*sin(q3) + 8.0*dq1*dq2*cos(q2)^2*cos(q3)^2 + 8.0*dq1*dq3*cos(q2)^2*cos(q3)^2 + 4.0*ddq1*cos(q2)*cos(q3)^2*sin(q2) + 4.0*ddq1*cos(q2)^2*cos(q3)*sin(q3) - 8.0*dq1*dq2*cos(q2)*cos(q3)*sin(q2)*sin(q3) - 8.0*dq1*dq3*cos(q2)*cos(q3)*sin(q2)*sin(q3), dq2^2*sin(q2)*sin(q3) - 1.0*ddq2*cos(q3)*sin(q2) - 1.0*ddq3*cos(q2)*sin(q3) - 1.0*ddq3*cos(q3)*sin(q2) - 1.0*dq2^2*cos(q2)*cos(q3) - 1.0*dq3^2*cos(q2)*cos(q3) - 1.0*ddq2*cos(q2)*sin(q3) + dq3^2*sin(q2)*sin(q3) - 2.0*dq2*dq3*cos(q2)*cos(q3) + 2.0*dq2*dq3*sin(q2)*sin(q3), 0.73*ddq1 - 0.558*ddq1*sin(q3) - 0.73*ddq1*cos(q2)^2 - 0.73*ddq1*cos(q3)^2 + 1.46*ddq1*cos(q2)^2*cos(q3)^2 + 0.558*ddq1*cos(q2)^2*sin(q3) - 0.558*dq1*dq2*cos(q3) - 0.558*dq1*dq3*cos(q3) + 1.116*dq1*dq2*cos(q2)^2*cos(q3) + 0.558*dq1*dq3*cos(q2)^2*cos(q3) + 0.558*ddq1*cos(q2)*cos(q3)*sin(q2) + 1.46*dq1*dq2*cos(q2)*sin(q2) + 1.46*dq1*dq3*cos(q2)*sin(q2) + 1.46*dq1*dq2*cos(q3)*sin(q3) + 1.46*dq1*dq3*cos(q3)*sin(q3) - 2.92*dq1*dq2*cos(q2)*cos(q3)^2*sin(q2) - 2.92*dq1*dq2*cos(q2)^2*cos(q3)*sin(q3) - 2.92*dq1*dq3*cos(q2)*cos(q3)^2*sin(q2) - 2.92*dq1*dq3*cos(q2)^2*cos(q3)*sin(q3) - 1.46*ddq1*cos(q2)*cos(q3)*sin(q2)*sin(q3) - 1.116*dq1*dq2*cos(q2)*sin(q2)*sin(q3) - 0.558*dq1*dq3*cos(q2)*sin(q2)*sin(q3), 2.0*dq1*dq2 - 4.0*dq1*dq2*cos(q2)^2 - 2.0*ddq1*cos(q2)*sin(q2), sin(q2)*dq2^2 - 1.0*ddq2*cos(q2), ddq1 - 1.0*ddq1*cos(q2)^2 - 1.0*ddq1*cos(q3)^2 + 2.0*ddq1*cos(q2)^2*cos(q3)^2 + 2.0*dq1*dq2*cos(q2)*sin(q2) + 2.0*dq1*dq3*cos(q2)*sin(q2) + 2.0*dq1*dq2*cos(q3)*sin(q3) + 2.0*dq1*dq3*cos(q3)*sin(q3) - 4.0*dq1*dq2*cos(q2)*cos(q3)^2*sin(q2) - 4.0*dq1*dq2*cos(q2)^2*cos(q3)*sin(q3) - 4.0*dq1*dq3*cos(q2)*cos(q3)^2*sin(q2) - 4.0*dq1*dq3*cos(q2)^2*cos(q3)*sin(q3) - 2.0*ddq1*cos(q2)*cos(q3)*sin(q2)*sin(q3), cos(q2)*dq2^2 + ddq2*sin(q2), ddq1 - 1.0*ddq1*cos(q2)^2 + 2.0*dq1*dq2*cos(q2)*sin(q2),           0, 0.558*dq1*dq2*sin(q3) - 0.73*dq1*dq3 - 0.558*ddq1*cos(q3) - 0.73*dq1*dq2 + 0.558*dq1*dq3*sin(q3) + 1.46*dq1*dq2*cos(q2)^2 + 1.46*dq1*dq2*cos(q3)^2 + 1.46*dq1*dq3*cos(q2)^2 + 1.46*dq1*dq3*cos(q3)^2 + 0.73*ddq1*cos(q2)*sin(q2) + 0.73*ddq1*cos(q3)*sin(q3) + 0.558*ddq1*cos(q2)^2*cos(q3) - 1.116*dq1*dq2*cos(q2)^2*sin(q3) - 0.558*dq1*dq3*cos(q2)^2*sin(q3) - 0.558*ddq1*cos(q2)*sin(q2)*sin(q3) - 2.92*dq1*dq2*cos(q2)^2*cos(q3)^2 - 2.92*dq1*dq3*cos(q2)^2*cos(q3)^2 - 1.46*ddq1*cos(q2)*cos(q3)^2*sin(q2) - 1.46*ddq1*cos(q2)^2*cos(q3)*sin(q3) - 1.116*dq1*dq2*cos(q2)*cos(q3)*sin(q2) - 0.558*dq1*dq3*cos(q2)*cos(q3)*sin(q2) + 2.92*dq1*dq2*cos(q2)*cos(q3)*sin(q2)*sin(q3) + 2.92*dq1*dq3*cos(q2)*cos(q3)*sin(q2)*sin(q3),   0,   0, dq1,    0, ddq1*cos(q2)^2 - 2.0*dq1*dq2*cos(q2)*sin(q2),   0,   0, 1.0, 1.116*dq1*dq2*cos(q2)^2 - 0.558*dq1*dq2 + 0.558*ddq1*cos(q2)*sin(q2), 0.558*ddq1 - 0.558*ddq1*cos(q2)^2 + 1.116*dq1*dq2*cos(q2)*sin(q2),         0, sign(dq1),         0];...
% [                                                                                                                                                                                                               ddq1*sin(q2)*sin(q3) - 1.0*ddq1*cos(q2)*cos(q3),                                                                                                                                                                                                                                                                                         2.0*dq1^2*cos(q2)^2 + 2.0*dq1^2*cos(q3)^2 - 1.0*dq1^2 - 4.0*dq1^2*cos(q2)^2*cos(q3)^2 + 4.0*dq1^2*cos(q2)*cos(q3)*sin(q2)*sin(q3),                                                                                                                                                                                                                         - 1.0*ddq1*cos(q2)*sin(q3) - 1.0*ddq1*cos(q3)*sin(q2),                                                                                                                                                                                                                                                                                                                              0.73*ddq2 + 0.73*ddq3 + 0.981*cos(q2)*cos(q3) - 0.981*sin(q2)*sin(q3) - 0.558*ddq2*sin(q3) - 0.279*ddq3*sin(q3) + 0.279*dq1^2*cos(q3) - 0.279*dq3^2*cos(q3) - 0.558*dq1^2*cos(q2)^2*cos(q3) - 0.73*dq1^2*cos(q2)*sin(q2) - 0.73*dq1^2*cos(q3)*sin(q3) - 0.558*dq2*dq3*cos(q3) + 1.46*dq1^2*cos(q2)*cos(q3)^2*sin(q2) + 1.46*dq1^2*cos(q2)^2*cos(q3)*sin(q3) + 0.558*dq1^2*cos(q2)*sin(q2)*sin(q3),                                2.0*dq1^2*cos(q2)^2 - 1.0*dq1^2,                -1.0*ddq1*cos(q2),                                                                                                                                                                                                                                                                                2.0*dq1^2*cos(q2)*cos(q3)^2*sin(q2) - 1.0*dq1^2*cos(q3)*sin(q3) - 1.0*dq1^2*cos(q2)*sin(q2) + 2.0*dq1^2*cos(q2)^2*cos(q3)*sin(q3),                 ddq1*sin(q2),                              -1.0*dq1^2*cos(q2)*sin(q2), ddq2 + ddq3,                                                                                                                                                                                                                                                                                                                                                          0.279*dq3^2*sin(q3) - 0.279*dq1^2*sin(q3) - 0.981*cos(q2)*sin(q3) - 0.981*cos(q3)*sin(q2) - 0.73*dq1^2*cos(q2)^2 - 0.73*dq1^2*cos(q3)^2 - 0.558*ddq2*cos(q3) - 0.279*ddq3*cos(q3) + 0.365*dq1^2 + 0.558*dq2*dq3*sin(q3) + 0.558*dq1^2*cos(q2)^2*sin(q3) + 1.46*dq1^2*cos(q2)^2*cos(q3)^2 + 0.558*dq1^2*cos(q2)*cos(q3)*sin(q2) - 1.46*dq1^2*cos(q2)*cos(q3)*sin(q2)*sin(q3), dq2,   0,   0, ddq2,                        dq1^2*cos(q2)*sin(q2), 1.0,   0,   0,                  0.981*cos(q2) - 0.558*dq1^2*cos(q2)^2 + 0.279*dq1^2,        - 0.558*cos(q2)*sin(q2)*dq1^2 + 0.558*ddq2 + 0.981*sin(q2),         0,         0, sign(dq2)];...
% [                                                                                                                                                                                                               ddq1*sin(q2)*sin(q3) - 1.0*ddq1*cos(q2)*cos(q3),                                                                                                                                                                                                                                                                                         2.0*dq1^2*cos(q2)^2 + 2.0*dq1^2*cos(q3)^2 - 1.0*dq1^2 - 4.0*dq1^2*cos(q2)^2*cos(q3)^2 + 4.0*dq1^2*cos(q2)*cos(q3)*sin(q2)*sin(q3),                                                                                                                                                                                                                         - 1.0*ddq1*cos(q2)*sin(q3) - 1.0*ddq1*cos(q3)*sin(q2),                                                                                                                                                                                                                                                                                                                                                                           0.73*ddq2 + 0.73*ddq3 + 0.981*cos(q2)*cos(q3) - 0.981*sin(q2)*sin(q3) - 0.279*ddq2*sin(q3) + 0.279*dq1^2*cos(q3) + 0.279*dq2^2*cos(q3) - 0.279*dq1^2*cos(q2)^2*cos(q3) - 0.73*dq1^2*cos(q2)*sin(q2) - 0.73*dq1^2*cos(q3)*sin(q3) + 1.46*dq1^2*cos(q2)*cos(q3)^2*sin(q2) + 1.46*dq1^2*cos(q2)^2*cos(q3)*sin(q3) + 0.279*dq1^2*cos(q2)*sin(q2)*sin(q3),                                                              0,                                0,                                                                                                                                                                                                                                                                                2.0*dq1^2*cos(q2)*cos(q3)^2*sin(q2) - 1.0*dq1^2*cos(q3)*sin(q3) - 1.0*dq1^2*cos(q2)*sin(q2) + 2.0*dq1^2*cos(q2)^2*cos(q3)*sin(q3),                            0,                                                       0, ddq2 + ddq3,                                                                                                                                                                                                                                                                                                                                                                                                       0.365*dq1^2 - 0.279*dq2^2*sin(q3) - 0.981*cos(q2)*sin(q3) - 0.981*cos(q3)*sin(q2) - 0.73*dq1^2*cos(q2)^2 - 0.73*dq1^2*cos(q3)^2 - 0.279*ddq2*cos(q3) - 0.279*dq1^2*sin(q3) + 0.279*dq1^2*cos(q2)^2*sin(q3) + 1.46*dq1^2*cos(q2)^2*cos(q3)^2 + 0.279*dq1^2*cos(q2)*cos(q3)*sin(q2) - 1.46*dq1^2*cos(q2)*cos(q3)*sin(q2)*sin(q3),   0, dq3,   0,    0,                                            0,   0, 1.0,   0,                                                                    0,                                                                 0, sign(dq3),         0,         0]];
%param = sym("param", [size(h_b,2),1]);
param = XB1_ols;
TAU = expand(A*param);
G = subs(TAU, {dq1 dq2 dq3, ddq1 ddq2 ddq3}, {0, 0, 0, 0, 0, 0});
CF = subs(TAU, {ddq1 ddq2 ddq3}, {0, 0, 0}) - G;
M_all = TAU - G - CF;
M = equationsToMatrix(M_all, {ddq1, ddq2, ddq3});
% ccode(TAU,'File','ccode/tau.c','Comments','Version: 1.0')
% ccode(G,'File','ccode/g.c','Comments','Version: 1.0')
% ccode(CF,'File','ccode/cf.c','Comments','Version: 1.0')
% ccode(M,'File','ccode/m.c','Comments','Version: 1.0')

%subs(M, {param(1), param(2), param(3), param(4), param5, param6, param7, param8, param9, param10, param11, param12, param13, param14, param15, param16, param17, param18, param19, param20, param21, param22, param23, param24},...
%     {XB1_wls(1), XB1_wls(2), XB1_wls(3), XB1_wls(4), XB1_wls(5), XB1_wls(6), XB1_wls(7), XB1_wls(8), XB1_wls(9), XB1_wls(10), XB1_wls(11), XB1_wls(12), XB1_wls(13), XB1_wls(14), XB1_wls(15), XB1_wls(16), XB1_wls(17), XB1_wls(18), XB1_wls(19), XB1_wls(20), XB1_wls(21), XB1_wls(22), XB1_wls(23), XB1_wls(24)})

%%
%A_t = subs(A, {ddq1, ddq2, ddq3, dq1, dq2, dq3, q2, q3}, {ddq1t, ddq2t, ddq3t, dq1t, dq2t, dq3t, q2t, q3t});

A_t = subs(A, {ddq1, ddq2, ddq3, dq1, dq2, dq3, q2, q3}, {diff(q1t(t), t, t),...
    diff(q2t(t), t, t), diff(q3t(t), t, t), diff(q1t(t),t), diff(q2t(t),t), diff(q3t(t),t),...
    q2t, q3t});
TAU_t = expand(A_t*param);
G_t = subs(TAU_t, {dq1t dq2t dq3t, ddq1t ddq2t ddq3t}, {0 0 0, 0 0 0});
CF_t = subs(TAU_t, {ddq1t ddq2t ddq3t}, {0 0 0});
M_all_t = TAU_t - G_t - CF_t;

%% Make ODE function

TAU_t = simplify(TAU_t);
[eqn,vars,Rnew] = reduceDifferentialOrder(TAU_t == 0, [q1t q2t q3t]);
[Mass,F] = massMatrixForm(eqn,vars);

M = simplify(Mass);
F = simplify(F);

f = M\F; % generate state-space form 
