clear variables
close all
%% ECE6320 Homework 2, Kade Howes

% 2.3
g = 9.8; m = 1/9.8; l = 0.25; b = 1;
% a
Aa = [0 1; g/l (-b)/(m*l^2)];
Ba = [0; 1/(m*l^2)];
[va, ea] = eig(Aa);
eacheck = roots([1 156.8 -39.2]);
% b
g = 9.8; m = 1/9.8; l = 0.25; b = 1;
Ab = [0 1; -g/l (-b)/(m*l^2)];
Bb = [0; 1/(m*l^2)];
[vb, eb] = eig(Ab);
ebcheck = roots([1 156.8 39.2]);
% b
g = 9.8; m = 1/9.8; l = 0.25; b = 1;
Ac = [0 1; g/(l*sqrt(2)) (-b)/(m*l^2)];
Bc = [0; 1/(m*l^2)];
[vc, ec] = eig(Ac);
eccheck = roots([1 156.8 -Ac(2,1)]);


%% Simulation of Linearized System
ueq = @(t,x) -g*l*m/(sqrt(2)); % equilibrium control input
xeq = [pi/4; 0];
SimulationScriptLinear(Ac, Bc, vc(:,2), ueq, xeq, 'Linear: Eigenvector - Negative Real Eigenvalue');
SimulationScriptLinear(Ac, Bc, vc(:,1), ueq, xeq, 'Linear: Eigenvector - Positive Real Eigenvalue');
SimulationScriptLinear(Ac, Bc, [(pi/4-0.05);0], ueq, xeq, 'Linear: theta = pi/4-0.05, thetadot = 0');

SimulationScriptLinear(Ac, Bc, [0;0], @(t,x) -g*l*m/(sqrt(2)), xeq, 'Linear: theta = pi/4, thetadot = 0');

%% Simulation of Nonlinear System
f = @(t,x,u) [x(2); (g*sin(x(1)/l)-(b*x(2)/(m*l^2))+(u/(m*l^2)))]
SimulationScriptNonlinear(f, vc(:,2)+xeq, ueq, 'Nonlinear: Eigenvector - Negative Real Eigenvalue');
SimulationScriptNonlinear(f, vc(:,1)+xeq, ueq, 'Nonlinear: Eigenvector - Positive Real Eigenvalue');
SimulationScriptNonlinear(f, [(pi/4-0.05);0]+xeq, ueq, 'Nonlinear: theta = pi/4-0.05, thetadot = 0');

SimulationScriptNonlinear(f, [pi/4;0]+xeq, ueq, 'Nonlinear: theta = pi/4, thetadot = 0');



