clear;
close all;
% Kade Howes, ECE6320, HW12

%% 2.1
Modified2_1;

%% LQR (Bryson's Method)
S = diag([1/10^2 1/1^2 1/2^2]);
LQRBryson(S, 'first');
S = diag([10 1 0.25]);
LQRBryson(S, 'second');

%% 2.3
A = [0 1; 1 0]; B = [0; 1];
Q = [1 0; 0 0]; R = [1];
syms p11 p12 p22;
P = [p11 p12; p12 p22];
eqn = A'*P + P*A + Q - P*B*inv(R)*B'*P;
p12_sol = 1+sqrt(2);
p22_sol = sqrt(2+2*sqrt(2));
p11_sol = 2*sqrt(1+sqrt(2));

P_sol = [p11_sol p12_sol; p12_sol p22_sol];
check_ricatti = double(subs(eqn,P,P_sol));
K_sol = [1+sqrt(2) sqrt(2+2*sqrt(2))];
K_check = lqr(A,B,Q,R);

%% 2.4
Ex2_4;

%% 2.5
A = [0 1; 0 0]; B = [0; 1];
Q = [4 0; 0 0]; R = [1];
K = lqr(A,B,Q,R);
eigClose = eig(A-B*K);

%% Output Feedback Control Design using Decompositions
decomposition_problem;




