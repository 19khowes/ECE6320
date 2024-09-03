function S05_L04_Tracking_Simple()
close all;

% Setup system dynamics (simple spring mass system)
A = [0 1; -4 -1];
B = [0; 10];
E = [1;0];
C = [1 0];

% Set simulation time variables
t0 = 0;
dt = 0.1;
tf = 10;
t_vec = [t0:dt:tf];

% Calculate feedback matrix using LQR
Q = diag([1, 0]);
R = 1;
K = lqr(A, B, Q, R);

% Calculate reference trajectory
ar  = [-6; 0];
r0 = [1;  3];
rSol = ode45(@(t,x)referenceDynamics(t, x, ar), [t0:dt:tf], r0);

% Calculate disturbance
aw  = [0; 0];
w0 = [0; 0]; % No disturbance
w0 = [1; 0]; % Constant disturbance
wSol = ode45(@(t,x)referenceDynamics(t, x, aw), [t0:dt:tf], w0);

% Simulate state
x0 = [0;0];
[t_mat, x_mat] = ode45(@(t, x)stateDynamics(t, x, K, rSol, wSol), [t0:dt:tf], x0 );
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
    function xdot = stateDynamics(t, x, K, rSol, wSol)
        % Calculate reference and disturbance
        r = deval(rSol, t);
        w = deval(wSol, t);
                
        % Get reference derivatives
        rdot = referenceDynamics(t, r, ar);
        rddot = rdot(end);
        
        % Calculate control
        z = x - r;
        u = -K*z + .1*(4*r(1) + r(2) + rddot);
        
        % Calculate state dynamics
        xdot = A*x + B*u + E*w(1);
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