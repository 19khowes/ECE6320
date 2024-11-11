clear;
close all;
% Kade Howes, ECE6320, HW11

% 16.2
A = [1 0 0; 1 1 0; -2 1 1];
b = [2; 0; 0];
c = [1 0 1];
Omega = obsv(A,c);
rOmega = rank(Omega);
syms k1 k2 k3 lam
k = [k1; k2; k3];
chi = charpoly(A+k*c);
phi = (lam+1)^3;
coef = coeffs(phi,lam); % get desired coefficients
ksol = solve(chi==coef); % use desired to get k's
newk = [ksol.k1; ksol.k2; ksol.k3];
newchi = subs(chi,k,newk); % check
eigAkc = eig(A+newk*c); % check


% Output feedback control design
A = [-5 0 0; 0 3 0; 0 2 1];
B = [0 0; 0 1; 1 0];
C = [0 0 1; -1 0 0];
xd = [0;1;0];
syms u1 u2;
uff = [u1;u2];
eqn = B*uff + A*xd == 0;
uff_sol = solve(eqn,uff);
uff_sol = double([uff_sol.u1, uff_sol.u2]);
rGamma = rank(ctrb(A,B)); % check controllability
maxX = [5; 3; 1];
maxU = [10; 1];
Q = diag(1./(maxX.^2)); % Bryson's method
R = diag(1./(maxU.^2));
K = lqr(A,B,Q,R)
L = place(A',C',[-100, -200, -300])'



