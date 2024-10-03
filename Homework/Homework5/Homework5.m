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

%% extra credit
close all;
% Linearized Stability
A1 = [0 1; -39.2 -156.8];
Q = eye(height(A1));
P = lyap(A1', Q);
eigP = eig(P);
mu = -1*(min(eig(Q))/max(eigP));

% finding a better Q (brute force style check for fun, didn't work?)
% bestMu = mu;
% bestQ = Q;
% bestP = P;
% for i = 1:1:100
%     for j = 1:1:100
%         for k = 1:1:100
%             Q = [i j; 0 k];
%             P = lyap(A1', Q);
%             new_mu = -1*min(eig(Q))/max(eig(P));
%             if (new_mu < bestMu && (min(eig(P)) > 0) && (sum((real(eig(P))~=eig(P))) == 0) )
%                 bestMu = new_mu;
%                 bestQ = Q;
%                 bestP = P;
%             end
%         end
%     end
% end
% mu = bestMu;
% P = bestP;
% Q = bestQ;

% try a few Q's
Q = [1 1; 0 1];
P = lyap(A1', Q);
mu = -1*min(eig(Q))/max(eig(P));
eig(P)

% b Simulate nonlinear system
xeq = [pi; 0];
dx_0 = [0.01; 0.75];
% convert dx_0 to x_0 for simulation
x_0 = dx_0 + xeq;
[tvec, xvec] = ode45(@(t,x) fnonlinear(t,x), [0 20], x_0);
xvec = xvec';

% back to dx
xvec = xvec - xeq;

v_t0 = dx_0'*P*dx_0; % change to x_0?
bounding = 1/(min(eig(P))) .* exp(mu.*tvec) .* v_t0;
dx_norm = zeros(1, length(tvec));
for k = 1:length(tvec)
    dx_norm(k) = norm(xvec(:,k))^2;
end

figure
plot(tvec, dx_norm, 'LineWidth', 3); hold on;
plot(tvec, bounding, 'LineWidth', 3);
title("Nonlinear System");
xlabel("time (t) [s]");
ylabel("||x||^2");
legend(["Nonlinear dynamics", "Bounding function"]);
% c Simulate linear system
[tvec, xvec] = ode45(@(t,x) flinear(t,x,A1), [0 20], dx_0);
xvec = xvec';

v_t0 = dx_0'*P*dx_0;
bounding = 1/(min(eig(P))) .* exp(mu.*tvec) .* v_t0;
dx_norm = zeros(1, length(tvec));
for k = 1:length(tvec)
    dx_norm(k) = norm(xvec(:,k))^2;
end

figure
plot(tvec, dx_norm, 'LineWidth', 3); hold on;
plot(tvec, bounding, 'LineWidth', 3);
title("Linear System");
xlabel("time (t) [s]");
ylabel("||x||^2");
legend(["Linear dynamics", "Bounding function"]);


%% Lyapunov equation
% Problem 1
% check if Q is positive definite
Q_1 = [1 2; 0 3];
eigQ_1 = eig(Q_1);
% check if resulting Ps are positive definite
P1 = [4.5 -0.866; -0.866 3.5]; P2 = [4.5 -0.866; 0.866 3.5];
eigP1_1 = eig(P1);
eigP2_1 = eig(P2);
% compute mu
mu_1 = -1*(min(eig(Q_1))/max(eig(P1)))
% Problem 2
% check if Q is positive definite
Q_2 = [1 1 ; 2 3];
eigQ_2 = eig(Q_2)
% check if resulting Ps are positive definite
P1 = [5 4; 3 2]; P2 = [-1 7; 0 3];
eigP1_2 = eig(P1)
eigP2_2 = eig(P2)


function xdot = flinear(t, x, A)
    xdot = A*x;
end

function xdot = fnonlinear(t,x)
    g = 9.8; m = 1/9.8; l = 0.25; b = 1;
    xdot = [x(2); g/l*sin(x(1))-b/(m*l^2)*x(2)];
end