function [t,y] = vdpDynamics(y0,T)
    options = odeset('RelTol',1e-6);
    
    [t,y] = ode45(@vdp1,[0 T],y0,options);

end
