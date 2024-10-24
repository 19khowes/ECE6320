clear;
close all;
%% Kade Howes, ECE6320, HW8

% 12.6
syms w real;
A = [0 1 0 0; 3*w^2 0 0 2*w; 0 0 0 1; 0 -2*w 0 1];
B = [0 0; 1 0; 0 0; 0 1];
Gamma = [B A*B A^2*B A^3*B];
rGamma_12_6 = rank(Gamma);

% 12.7
syms x1 x2 x3 u1 u2;
f = [-x1+u1; -x2+u2; x2*u1-x1*u2];
x = [x1; x2; x3];
u = [u1; u2];
dfdx = jacobian(f,x);
dfdu = jacobian(f,u);
A1 = subs(dfdx, [x1 x2 x3 u1 u2], [0 0 0 0 0]);
B1 = subs(dfdu, [x1 x2 x3 u1 u2], [0 0 0 0 0]);
Gamma = [B1 A1*B1 A1^2*B1];
rGamma_12_7_a = rank(Gamma);
A2 = subs(dfdx, [x1 x2 x3 u1 u2], [1 1 1 1 1]);
B2 = subs(dfdu, [x1 x2 x3 u1 u2], [1 1 1 1 1]);
Gamma = [B2 A2*B2 A2^2*B2];
rGamma_12_7_b = rank(Gamma);

% 13.1 
A = [-1 0; 0 -1]; B = [-1; 1]; C = [1 0; 0 1]; D = [2; 1];
Gamma = [B A*B];
rGamma_13_1 = rank(Gamma);
v = orth(Gamma);
w = null(Gamma');
T = [v w];
Ahat = T^-1*A*T;
Bhat = T^-1*B;
Chat = C*T;
Dhat = D;
syms s;
Ahat = round(Ahat);
G1 = simplify(C*((s*eye(2)-A)^-1)*B+D);
G2 = simplify(Chat*((s*eye(2)-Ahat)^-1)*Bhat+Dhat);


% Control Design Review
% A
A = [0 -2; 2 0]; B = [0; 1];
Gamma = [B A*B];
rGamma_A = rank(Gamma);
syms k1 k2;
Abar = A-B*[k1 k2];
Abar = subs(Abar,[k1 k2],[0 4]);
eigAbar = eig(Abar);
% B
A = [1 0; 0 -1]; B = [0; 1];
Gamma = [B A*B];
rGamma_B = rank(Gamma);
v = orth(Gamma);
w = null(Gamma');
T = [v w];
Ahat = T^-1*A*T;
Bhat = T^-1*B;
% C
A = [-1 0; 0 1]; B = [0; 1];
Gamma = [B A*B];
rGamma_C = rank(Gamma);
v = orth(Gamma);
w = null(Gamma');
T = [v w];
Ahat = T^-1*A*T;
Bhat = T^-1*B;

syms k1 k2;
Abar = A-B*[k1 k2];
Abar = subs(Abar,[k1 k2],[0 4]);
eigAbar = eig(Abar);

%% extra credit
syms a1 a2 a3 A B
C = [B A*B A^2*B];
T1 = [1 a1 a2; 0 1 a1; 0 0 1];
lh = A*C*T1;
rh1 = [-a1 -a2 -a3; 1 0 0; 0 1 0];
rh = C*T1*rh1;

A = [6 4 1; -5 -4 0; -4 -3 -1];
B = [1; -1; -1];
chiA = poly(A);
a1 = chiA(2);
a2 = chiA(3);
C = ctrb(A,B);
T = C*[1 a1 a2; 0 1 a1; 0 0 1];
Ahat = T^-1*A*T;
Bhat = T^-1*B;




