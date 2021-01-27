%%% Neron Oscillator Simulation

%% Van del Pol
y0= [2;0];
T = 20;
[t,y] = vdpDynamics(y0,T);
% plot
plot(t,y(:,1),'-o',t,y(:,2),'-o')
title('Solution of van der Pol Equation (\mu = 1) with ODE45');
xlabel('Time t');
ylabel('Solution y');
legend('y_1','y_2');

%% Hindmarsh-Rose
X0= [1;1;1];
T = 20;
I = 10;
param.a = 1; param.b = 3; param.c = 1; param.d = 5; 
param.r = 1e-3; param.s = 4; param.V0 = -1.6;
[t,X] = HRDynamics(X0,T,I,param);
% plot
plot(t,X(:,1),'-o',t,X(:,2),'-o',t,X(:,3),'-o')
title('Solution of Hindmarsh-Rose Equation with ODE45');
xlabel('Time t');
ylabel('Solution y');
legend('X_1','X_2','X_3');

%% Hodgkin-Huxley
X0= [1;1;1;1];
T = 20;
I = -50;
param.VNa = 50; param.VK = -77; param.VL = -54.4; param.gNa = 120; 
param.gK = 36; param.gL = 0.3; param.Ib = 10; param.c = 1;
[t,X] = HHDynamics(X0,T,I,param);
% plot
plot(t,X(:,1),'-o',t,X(:,2),'-o',t,X(:,3),'-o',t,X(:,3),'-o')
title('Solution of Hodgkin-Huxley Equation with ODE45');
xlabel('Time t');
ylabel('Solution y');
legend('X_1','X_2','X_3','X_4');