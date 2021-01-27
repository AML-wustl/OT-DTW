function [dX] = HHFunction(t,X,I,param)
    %parameters
    if(nargin < 4)
    VNa = 50; VK = -77; VL = -54.4; gNa = 120; gK = 36; gL = 0.3; Ib = 10; c = 1;
    elseif(nargin == 4)
        VNa = param.VNa; VK = param.VK; VL = param.VL; gNa = param.gNa;
        gK = param.gK; gL = param.gL; Ib = param.Ib; c = param.c;
    end
    
    if(nargin < 3) 
        I = 0; 
    end
    V = X(1); m = X(2); h = X(3); n = X(4); 
    
    % HH model equations
    am = 0.1*(V + 40)./(1 - exp(-(V + 40)/10));
    bm = 4*exp(-(V + 65)/18);
    ah = 0.07*exp(-(V + 65)/20);
    bh = 1./(1 + exp(-(V + 35)/10));
    an = 0.01*(V + 55)./(1 - exp(-(V + 55)/10));
    bn = 0.125*exp(-(V + 65)/80);
    dV = (Ib + I - gNa*h.*(V - VNa).*(m.^3) - gK*(V - VK).*(n.^4) - gL*(V - VL))/c;
    dm = am.*(1-m) - bm.*m;
    dh = ah.*(1-h) - bh.*h;
    dn = an.*(1-n) - bn.*n;
    
    dX = [dV; dm; dh; dn];
end