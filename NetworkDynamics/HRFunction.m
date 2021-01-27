function [dX] = HRFunction(t,X,I,param)
    if(nargin < 4)
        a = 1; b = 3; c = 1; d = 5; r = 1e-3; s = 4; V0 = -1.6; 
    elseif(nargin == 4)
        a = param.a; b = param.b; c = param.c; d = param.d; r = param.r; 
        s = param.s; V0 = param.V0; 
    end
    if(nargin < 3) 
        I = 0; 
    end
    
    V = X(1); 
    n = X(2); 
    h = X(3); 
    %HR model equation
    dV = n - a*V^3 +b*V^2 - h + I;
    dn = c - d*V^2 - n;
    dh = r*(s*(V - V0) - h);
    
    dX = [dV; dn; dh];
end