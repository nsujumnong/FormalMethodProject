function model = makeODE(eqn)
% This function is based on the ODE function generation from Radian's
% impedance control code

%     dyn_eqn=ctrl.u-eqn;

%     dyn_eqn=simplify(dyn_eqn);
    %dyn_eqn = eqn==u;

    [eqs,vars] = reduceDifferentialOrder(eqn==0,[q1t,q2t,q3t]);
    [Mass,F] = massMatrixForm(eqs,vars);

    % 
    % F = vpa(F,2);
    % M = vpa(M,2);


    F=simplify(F);
    M=simplify(Mass);

    model = M\F; too slow doesnt work

end