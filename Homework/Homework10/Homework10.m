clear;
close all;
% Kade Howes, ECE6320, HW10

% 15.4
% a
A = [-2 0 0 0; 0 -2 0 0; 1 0 0 0; 0 1 0 0];
B = [1 0; 0 1; 0 0; 0 0];
C = [1 1 2 0];
D = [1 0];
syms s;
G = simplify(C*(s*eye(height(A))-A)^-1 * B + D);
% b
rGamma = rank(ctrb(A,B));
rOmega = rank(obsv(A,C));

% 15.6
% a
A = [-1 10; 0 1]; B = [1; 0]; C = [1 1];
rGamma = rank(ctrb(A,B));
rOmega = rank(obsv(A,C));
% b & c
G = C*(s*eye(height(A))-A)^-1 * B;
% d
eigA = eig(A);

% 15.7
% a
syms x1 x2 u
x = [x1; x2];
f = [x2; 1-(u^2/x1^2)];
g = x1;
% b 
yeq = 1;
xeq = [yeq; 0];
f_sol = subs(f,x,xeq);
ueq = solve(f_sol,u);
% c
dfdx = jacobian(f,x);
dfdu = jacobian(f,u);
dgdx = jacobian(g,x);
A1 = subs(dfdx,[x; u],[xeq; ueq(2)]);
B1 = subs(dfdu,[x; u],[xeq; ueq(2)]);
C1 = subs(dgdx,[x; u],[xeq; ueq(2)]);
An1 = subs(dfdx,[x; u],[xeq; ueq(1)]);
Bn1 = subs(dfdu,[x; u],[xeq; ueq(1)]);
Cn1 = subs(dgdx,[x; u],[xeq; ueq(1)]);
% d
eigA1 = eig(A1);
eigAn1 = eig(An1);
% e
rGamma1 = rank(ctrb(A1,B1));
rGamman1 = rank(ctrb(An1,Bn1));
rOmega1 = rank(obsv(A1,C1));
rOmegan1 = rank(obsv(An1,Cn1));

% 17.4
% a
A = [2 1; 0 1]; B = [0; 1]; C = [1 1];
Gamma = ctrb(A,B);
rGamma = rank(Gamma);
Omega = obsv(A,C);
rOmega = rank(Omega);
% b
G = simplify(C*(s*eye(height(A))-A)^-1 * B);
% c
v = null(Omega);
u = orth(Omega');
T = [u v];
Ahat =T^-1*A*T;
Bhat = T^-1*B;
Chat = C*T;
A11hat = Ahat(1:rOmega,1:rOmega);
B1hat = Bhat(1:rOmega,:);
C1hat = Chat(:,1:rOmega);
Ghat = simplify(C1hat*(s*eye(height(A11hat))-A11hat)^-1 * B1hat);

% 17.4
% a
A = [1 1; 0 2]; B = [1; 1]; C = [1 0];
Gamma = ctrb(A,B);
rGamma = rank(Gamma);
Omega = obsv(A,C);
rOmega = rank(Omega);
% b
G = simplify(C*(s*eye(height(A))-A)^-1 * B);
% c
v = orth(Gamma);
u = null(Gamma');
T = [v u];
Ahat =T^-1*A*T;
Bhat = T^-1*B;
Chat = C*T;
A11hat = Ahat(1:rGamma,1:rGamma);
B1hat = Bhat(1:rGamma,:);
C1hat = Chat(:,1:rGamma);
Ghat = simplify(C1hat*(s*eye(height(A11hat))-A11hat)^-1 * B1hat);








