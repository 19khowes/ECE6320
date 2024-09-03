function S05_L04_Sinuousoid_Tracking()
close all;

% Setup system dynamics (simple spring mass system)
A = [0 1; -4 -1];
B = [0; 10];
E = [1;0];
C = [1 0];

% Set simulation time variables
t0 = 0;
dt = 0.1;
tf = 50;
t_vec = [t0:dt:tf];

% Create disturbance
omega_0 = 2;
aw  = [-omega_0^2; 0];
w0 = [1; 5 ];
wSol = ode45(@(t,x)referenceDynamics(t, x, aw), [t0:dt:tf], w0);

% Calculate reference trajectory
ar  = [-1; 0];
r0 = [1;  6];
rSol = ode45(@(t,x)referenceDynamics(t, x, ar), [t0:dt:tf], r0);

% Augmented matrix for calculating gain (note: ctrb has rank 4)
Abar = [0 1 0 0; omega_0^2 0 C; zeros(2,2) A];
Bbar = [0; 0; B];

% Calculate feedback matrix using LQR
Q = diag([1, 1, 1, 0]); % Note this is observable
R = 1;
K = lqr(Abar, Bbar, Q, R);

% Create the initial state and transformed state
x0 = [0; 0];
z0 = x0-r0;

% Create the initial aggregate state
phi_1_0 = 0;  % s^-2 e
phi_2_0 = 0;  % s^-1 e
phi_3_0 = zeros(2,1);  % s^-2 z
phi_4_0 = zeros(2,1); % s^-1 z
phi_5_0 = 0;  % s^-2 uhat
phi_6_0 = 0;  % s^-1 uhat
phi_7_0 = x0;  % state
phi0 = [phi_1_0; phi_2_0; phi_3_0; phi_4_0; phi_5_0; phi_6_0; phi_7_0];

% Create indices of variables
ind_phi_1 = 1;
ind_phi_2 = 2;
ind_phi_3 = [3:4];
ind_phi_4 = [5:6];
ind_phi_5 = 7;
ind_phi_6 = 8;
ind_phi_7 = [9:10];
ind_x = ind_phi_7;

% Simulate the state
[t_mat, phi_mat] = ode45(@(t, phi)stateDynamics(t, phi, K, rSol, wSol), [t0:dt:tf], phi0 );
t_mat = t_mat';
phi_mat = phi_mat';

% Extract data from phi_mat
x_mat = phi_mat(ind_x, :);

%% Plot results
% Get reference for plotting
err = zeros(length(t_mat));
r_mat = zeros(length(r0), length(t_mat));
for k = 1:length(t_mat)
    r_mat(:, k) = deval(rSol, t_mat(k));
    err(k) = x_mat(1,k) - r_mat(1,k);
end
ref = r_mat(1,:);

% Get disturbance for plotting
w_mat = zeros(length(r0), length(t_mat));
for k = 1:length(t_mat)
    w_mat(:, k) = deval(wSol, t_mat(k));
end
dist = w_mat(1,:);

% Plot states
figure;
subplot(2, 1, 1); % Plot x1
plot(t_vec, ref, 'r', 'linewidth', 3); hold on;
plot(t_vec, dist, 'm', 'linewidth', 2);
plot(t_mat, x_mat(1,:), 'b', 'linewidth', 2);
ylabel('y, r, and w');

subplot(2,1,2); % Plot error
plot(t_mat, err, 'b', 'linewidth', 2);
ylabel('Error (y - r)');
xlabel('Time(s)');


%% Dynamics function
    function phidot = stateDynamics(t, phi, K, rSol, wSol)
        % Extract states from phi
        phi_1 = phi(ind_phi_1);  % s^-2 e
        phi_2 = phi(ind_phi_2);  % s^-1 e
        phi_3 = phi(ind_phi_3);  % s^-2 z
        phi_4 = phi(ind_phi_4); % s^-1 z
        phi_5 = phi(ind_phi_5);  % s^-2 uhat
        phi_6 = phi(ind_phi_6);  % s^-1 uhat
        phi_7 = phi(ind_phi_7);  % state
        x = phi_7;
        
        % Calculate reference and disturbance
        r = deval(rSol, t);
        rdot = referenceDynamics(t, r, ar);
        rddot = rdot(end);
        rdot = rdot(1);
        w = deval(wSol, t);
        
        % Calculate z
        z = x - r;
                
        % Calculate control
        xi = [phi_1; phi_2; z + omega_0^2*phi_3];
        uhat = -K*xi - omega_0^2*phi_5;
        uff = .1*(4*r(1) + r(2) + rddot);
        u = uhat + uff;
        
        % Calculate aggregate state dynamics
        phi_1_dot = phi_2;
        phi_2_dot = C*z;
        phi_3_dot = phi_4;
        phi_4_dot = z;
        phi_5_dot = phi_6;
        phi_6_dot = uhat;
        phi_7_dot = A*x+B*u+E*w(1);
        phidot = [phi_1_dot; phi_2_dot; phi_3_dot; phi_4_dot; ...
                  phi_5_dot; phi_6_dot; phi_7_dot];
        
    end

end

function rdot = referenceDynamics(t, r, a)
    % This function implements a reference trajectory of the form:
    %       r^(p) = \sum_{i=0}^{p-1} a_i r^(i)
    %
    % Note that the book implements the reference trajectory as:
    %       r^(p) = \sum_{i=1}^p a_i r^(p-i)

    % Calculate the pth derivative
    r_p = a'*r;
    
    % Create vector
    rdot = zeros(length(a), 1);
    rdot(1:end-1) = r(2:end);
    rdot(end) = r_p;    
end