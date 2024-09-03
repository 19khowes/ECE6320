function S05_L04_Tracking_Integral()
close all;

% Setup system dynamics (simple spring mass system)
A = [0 1; -4 -1];
B = [0; 10];
E = [1;0];
C = [1 0];

% Augmented dynamics
Abar = [A zeros(2,1); C 0];
Bbar = [B; 0];
Ebar = [E; 0];

% Set simulation time variables
t0 = 0;
dt = 0.1;
tf = 10;
t_vec = [t0:dt:tf];

% Calculate feedback matrix using LQR
Q = diag([1, 0, 1]);
R = 1;
K = lqr(Abar, Bbar, Q, R);

% Calculate reference trajectory
ar  = [-2; 0];
r0 = [1;  3];
rSol = ode45(@(t,x)referenceDynamics(t, x, ar), [t0:dt:tf], r0);

% Calculate disturbance
aw  = [0; 0];
%w0 = [1; 0]; % Constant disturbance
w0 = [1; 1]; % Ramp disturbance
wSol = ode45(@(t,x)referenceDynamics(t, x, aw), [t0:dt:tf], w0);

% Simulate state
x0 = [0;0];
xbar0 = [x0; 0];
[t_mat, x_mat] = ode45(@(t, xbar)stateDynamics(t, xbar, K, rSol, wSol), [t0:dt:tf], xbar0 );
x_mat = x_mat';
t_mat = t_mat';

%% Plot results
% Get reference for plotting
err = zeros(1, length(t_mat));
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
    function xbardot = stateDynamics(t, xbar, K, rSol, wSol)
        % Calculate reference and disturbance
        r = deval(rSol, t);
        w = deval(wSol, t);
                
        % Get reference derivatives
        rdot = referenceDynamics(t, r, ar);
        rddot = rdot(end);
        
        % Calculate control
        x = xbar(1:2);  % Get original state
        z = x - r;  % Get shifted state
        xi = xbar - [r; 0]; % Get augmented state
        u = -K*xi + .1*(4*r(1) + r(2) + rddot);
        
        % Calculate state dynamics
        sig1_dot = C*z; % Get derivative of integral term
        xdot = A*x + B*u + E*w(1);
        xbardot = [xdot; sig1_dot];
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