%% Create the dynamics
% Symbolic variables for solving for dynamics
syms u real % Force applied to the cart
syms theta real % pendulum angle from vertical
syms thetadot real % time derivative of theta
syms thetaddot real % time derivative of thetadot
syms x real % cart position coordinate
syms xdot real % cart velocity
syms xddot real % cart acceleration

% Create xddot
t2 = cos(theta);
t3 = sin(theta);
xddot = ((u.*1.0e2-xdot.*1.0e1+t2.*t3.*1.47e2+t3.*thetadot.^2.*6.0).*(-1.0./5.0))./(t2.^2.*3.0-1.4e1);

% Create thetaddot
thetaddot = (t3.*3.43e2+t2.*u.*5.0e1-t2.*xdot.*5.0+t2.*t3.*thetadot.^2.*3.0)./(t2.^2.*3.0-1.4e1);

% Create dynamic function
f = [xdot; xddot; thetadot; thetaddot];

%% Create the partials
% Partial of dynamics with respect to the state
z = [x; xdot; theta; thetadot]; % z is the state
df_dz = jacobian(f, z);

% Partial of dynamics with respect to the input
df_du = jacobian(f,u);

%% Substitute in the state at the equilibrium
% set the equilibrium state
x = 0;
xdot = 0;
theta = pi;
thetadot = 0;

% evaluate the partials
disp('Partials calculated in matlab');
df_dz_eq = subs(df_dz)
df_du_eq = subs(df_du)

%% Evaluate the partials calculated by hand
% Create the partials by hand
df_dz_eq_hand = [0 1 0 0; 0 -10/55 147/55 0; 0 0 0 1; 0 -5/11 343/11 0];
df_du_eq_hand = [0; 100/55; 0; 50/11];

% Output the difference
disp('Difference between matlab and hand:');
err_df_dz = df_dz_eq - df_dz_eq_hand
err_df_du = df_du_eq - df_du_eq_hand
