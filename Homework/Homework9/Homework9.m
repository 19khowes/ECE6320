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
eigAbar = eig(Abar);
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
Abar = A-B*K;
eigAbar = eig(Abar);
S02_L03_PendulumEnergy([pi/4-0.1;0],xeq,ueq,K);
S02_L03_PendulumEnergy([pi-0.1;0],xeq,ueq,K);

% 2.8
% a
l = 1; m = 1; b = 0.1; g = 9.8;
A = [0 1; 0 -b/(m*l^2)];
B = [0; 1/(m*l^2)];
Gamma = ctrb(A,B);
rGamma = rank(Gamma);
K = lqr(A,B,diag(1:height(A)), 1); % control law for uhat
Abar = (A-B*K); % checking stabilization worked
eigAbar = eig(Abar);

% b

%% Cart-Pendulum
% These variables are described at http://ctms.engin.umich.edu/CTMS/index.php?example=InvertedPendulum&section=SystemModeling
M = 0.5; % Mass of cart
m = 0.2; % mass of pendulum
b = 0.1; % coefficient of friction for cart
l = 0.3; % length to pendulum center of mass
I = 0.006; % mass moment of inertia of the pendulum
g = 9.8; % Gravity constant

% Symbolic variables for solving for dynamics
syms u real % Force applied to the cart
syms theta real % pendulum angle from vertical
syms thetadot real % time derivative of theta
syms thetaddot real % time derivative of thetadot
syms x real % car position coordinate
syms xdot real % car velocity
syms xddot real % car acceleration
% Solve for the equations
eqn1 = (M+m)*xddot + b*xdot + m*l*thetaddot*cos(theta) - m*l*thetadot^2*sin(theta) - u;
eqn2 = (I+m*l^2)*thetaddot + m*g*l*sin(theta) + m*l*xddot*cos(theta);
soln = solve(eqn1, eqn2, xddot, thetaddot);
thetaddot = simplify(soln.thetaddot);
xddot = simplify(soln.xddot);

% state z
z = [x; xdot; theta; thetadot];

% define dynamics
f = [xdot; xddot; thetadot; thetaddot];

theta_sol = pi;
thetadot_sol = 0;
xdot_sol = 0;
x_sol = 0; % choose this, x_sol can be anything
f_sol = subs(f,theta,theta_sol);
f_sol = subs(f_sol,thetadot,thetadot_sol);
f_sol = subs(f_sol,xdot,xdot_sol);
u_sol = solve(f_sol,u);

% Create the jacobians
df_dz = jacobian(f, z); % Creates a jacobian by evaluating f wrt z
df_du = jacobian(f, u); % Create a jacobian with respect to u

% choose equilibrium
x = x_sol;
xdot = xdot_sol;
theta = theta_sol;
thetadot = thetadot_sol;
u = u_sol;

% Evaluate linearization at the equilibirum
df_dz_sol = simplify(subs(df_dz)); % The "subs" command will substitute all of the redefined variables into df_dz
df_du_sol = simplify(subs(df_du)); % Do a very similar thing to get the solution for df_du

% Get the A and B matrices (need to convert the symbolic matrices to double
% matrices)
A = double(subs(df_dz_sol)); % Takes the solution for df_dx and creates a matrix of doubles
B = double(subs(df_du_sol)); % Need to do the same thing for B

% evaluate stability
eigA = eig(A);

% create controller 
% evaluate controllability
Gamma = ctrb(A,B);
rGamma = rank(Gamma);

% create K
K = lqr(A,B,diag(1:rGamma),diag(1:width(B)));

% evaluate stability
Abar = A-B*K;
eigAbar = eig(Abar);

% plot 1
z0 = [0;0;pi-0.25;0];
zeq = @(t,x)[x_sol;xdot_sol;theta_sol;thetadot_sol];
ueq = @(t,x)u_sol;
PlotCart(z0,zeq,ueq,K)


% Symbolic variables for solving for dynamics
syms u real % Force applied to the cart
syms theta real % pendulum angle from vertical
syms thetadot real % time derivative of theta
syms thetaddot real % time derivative of thetadot
syms x real % car position coordinate
syms xdot real % car velocity
syms xddot real % car acceleration
syms t

% 2nd Equilibrium point
theta_sol = 3*pi/4;
thetadot_sol = 0;
% xdot_sol = 0;
% x_sol = 0; % choose this, x_sol can be anything
f_sol = subs(f,theta,theta_sol);
f_sol = subs(f_sol,thetadot,thetadot_sol);
% f_sol = subs(f_sol,xdot,xdot_sol);
u_sol = solve(f_sol(4)==0,u);
f_sol = subs(f_sol,u,u_sol);

% Create the jacobians
df_dz = jacobian(f, z); % Creates a jacobian by evaluating f wrt z
df_du = jacobian(f, u); % Create a jacobian with respect to u

% choose trajectory (solution)
x = (49/10)*t^2;
xdot = 49/5*t;
theta = theta_sol;
thetadot = thetadot_sol;
u = u_sol;

% Evaluate linearization at the equilibirum
df_dz_sol = simplify(subs(df_dz)); % The "subs" command will substitute all of the redefined variables into df_dz
df_du_sol = simplify(subs(df_du)); % Do a very similar thing to get the solution for df_du

% Get the A and B matrices (need to convert the symbolic matrices to double
% matrices)
A = double(subs(df_dz_sol)); % Takes the solution for df_dx and creates a matrix of doubles
B = double(subs(df_du_sol)); % Need to do the same thing for B

% evaluate stability
eigA = eig(A);

% create controller 
% evaluate controllability
Gamma = ctrb(A,B);
rGamma = rank(Gamma);

% create K
K = lqr(A,B,diag(1:rGamma),diag(1:width(B)));

% evaluate stability
Abar = A-B*K;
eigAbar = eig(Abar);

% plot
z0 = [0;0;pi-0.25;0];
zsol = @(t,x) [(49/10)*t^2;49/5*t;theta_sol;thetadot_sol];
usol = @(t,x) x(2)/10 + 343/50;
PlotCart(z0,zsol,usol,K)











