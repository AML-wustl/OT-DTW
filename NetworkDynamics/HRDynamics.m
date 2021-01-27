function [tt,X] = HRDynamics(X0,T,I,param)
    options = odeset('RelTol',1e-6);
    
    [tt,X] = ode45(@(t,X) HRFunction(t,X,I,param),[0 T],X0,options);

end
