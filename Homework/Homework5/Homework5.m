clear;
close all;

%% ECE 6320, Kade Howes, Homework 5

% 8.6
syms t real
syms t0 real
Phif = @(t) [exp(t)*cos(2*t) exp(-2*t)*sin(2*t); -exp(t)*sin(2*t) exp(-2*t)*cos(2*t)];
Phi_t = Phif(t);
Phi_t0 = Phif(t0);
Phi_t_t0 = simplify(Phi_t*(Phi_t0^-1));

% choosing tao = 0 to compute A
dPhi_t = simplify(diff(Phi_t,t));
A = simplify(dPhi_t*(Phi_t^-1));

% A = [(-1/2)+(3/2)*cos(4*t) 2-(3/2)*sin(4*t); -2-(3/2)*sin(4*t) (-1/2)-(3/2)*cos(4*t)]; % from homework answer, replace with own code
syms lam
Alam = lam*eye(height(A))-A;
charPolyA = det(Alam);
eigA = roots([1 1 2]);
eigAcheck = simplify(eig(A));

% 8.12
% a
syms x1 real
syms x2 real
x = [x1; x2];
f = [-x1+x1*(x1^2+x2^2); -x2+x2*(x1^2+x2^2)];
df_dx_a = jacobian(f,x);
% substitute equilibrium point of x
x1 = 0;
x2 = 0;
df_dx_aeq = subs(df_dx_a);
% b
syms x1 real
syms x2 real
syms alpha real
x = [x1; x2];
f = [x2; -1*alpha*x2-x1];
df_dx_b = simplify(jacobian(f,x));
% substitute equilibrium point of x
x1 = 0;
x2 = 0;
df_dx_beq = subs(df_dx_b);
% find eigenvalues
syms lam real
charPolyb = det(lam*eye(height(df_dx_beq))-df_dx_beq);

% 2.3
g = 9.8; m = 1/9.8; l = 0.25; b = 1;
% a
Aa = [0 1; g/l (-b/(m*l^2))];
eigAa = eig(Aa);
% b
Ab = [0 1; -g/l (-b/(m*l^2))];
eigAb = eig(Ab);
% c
Ac = [0 1; g/(l*sqrt(2)) (-b/(m*l^2))];
eigAc = eig(Ac);

% 2.4
k = 1;
A4c = [0 1; -k 0];
eigA4c = eig(A4c);
% eigenvalue check doesn't work, try pos def. Q
Q = eye(height(A4c));
% P = lyap(A4c',Q) % doesn't exist, fails

% 2.6
I = 1; b = 2;
A6b = [0 1; (g*m)/I (-b/I)];
eigA6b = eig(A6b);

% 2.7
A7d = [0 1 0; -1 0 0; 0 0 0];
eigA7d = eig(A7d);


