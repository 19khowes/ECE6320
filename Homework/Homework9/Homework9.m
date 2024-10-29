clear;
close all;
%% Kade Howes, ECE6320, HW9

% Control design using the controllability decomposition

% 1
control_decomposition();

% Section 2.6 - Control Design
% 2.3
g = 9.8; m = 1/9.8; l = 0.25; b = 1;
% a
xeq = [0;0]; ueq = 0;
A = [0 1; g/l (-b/(m*l^2))];
B = [0; 1/(m*l^2)];
xd = xeq;
xddot = [0;0];
% syms uff
% eqn = B*uff+A*xd-xddot==0;
% uff_sol = double(solve(eqn,uff));
% uff = 0; % TODO; replace with some solve for uff?
Gamma = ctrb(A,B);
rGamma = rank(Gamma);
Q = diag(1:rGamma);
R = diag(1:width(B));
K = lqr(A, B, Q, R); % control for controllable portion
Abar = A-B*K;
S02_L03_PendulumEnergy([0.1;0],xeq,ueq,K)
S02_L03_PendulumEnergy([(pi-0.1);0],xeq,ueq,K)

% c
xeq = [pi/4;0]; ueq = -g*l*m/sqrt(2);
A = [0 1; g/(l*sqrt(2)) (-b/(m*l^2))];
B = [0; 1/(m*l^2)];
xd = xeq;
xddot = [0;0];
syms uff
eqn = B*uff+A*xd-xddot==0;
uff_sol = double(solve(eqn,uff));
% uff = 0; % TODO; replace with some solve for uff?
Gamma = ctrb(A,B);
rGamma = rank(Gamma);
Q = diag(1:rGamma);
R = diag(1:width(B));
K = lqr(A, B, Q, R); % control for controllable portion
Abar = A-B*K
S02_L03_PendulumEnergy([pi/4-0.1;0],xeq,ueq,K)
S02_L03_PendulumEnergy([pi-0.1;0],xeq,ueq,K)












