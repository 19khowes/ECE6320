function S05_L06_L07_AircraftExample()
close all;
    % Seup for convergences to x = 0
    P = controlParamters(); % Create control
    u = @(t, x) -P.K*x; % Set the control handle
    x0 = rand(3,1); % Initial state 
    xhat0 = zeros(3,1);
    
    % Simulate the state forward in time the state
    dt = 0.001;
    t = [0:dt:1.5];
    [tmat, xmat] = ode45(@(t,x)dynamics(t,x,u(t, x(4:6)), P), t, [x0; xhat0]);
    tmat = tmat';
    xmat = xmat';
    
    % Calculate the energy
    len = length(tmat);
    E = zeros(1,len);
    u_mat = zeros(1,len);
    for k = 1:len
        u_mat(k) = u(tmat(k), xmat(4:6,k));
    end
    
    %% Plot the results
    fontsize = 12;
    
    % Plot the states
    figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    subplot(4,1,1);
    plot([tmat(1) tmat(end)], [P.xd(1) P.xd(1)], ':r', 'linewidth', 3); hold on;
    plot(tmat, xmat(4,:), 'g', 'linewidth', 2.5);
    plot(tmat, xmat(1,:), 'b', 'linewidth', 2);
    ylabel('$\theta(t)$', 'fontsize', fontsize, 'Interpreter','latex');
    set(gca, 'fontsize', fontsize);
    
    subplot(4,1,2);
    plot([tmat(1) tmat(end)], [P.xd(2) P.xd(2)], 'r:', 'linewidth', 3); hold on;
    plot(tmat, xmat(5,:), 'g', 'linewidth', 2.5);
    plot(tmat, xmat(2,:), 'b', 'linewidth', 2);
    ylabel('$\dot{\omega}(t)$', 'fontsize', fontsize, 'Interpreter','latex');
    set(gca, 'fontsize', fontsize);
    
    subplot(4,1,3);
    plot([tmat(1) tmat(end)], [P.xd(3) P.xd(3)], 'r:', 'linewidth', 3); hold on;
    plot(tmat, xmat(6,:), 'g', 'linewidth', 2.5);
    plot(tmat, xmat(3,:), 'b', 'linewidth', 2);
    ylabel('$\dot{\tau}(t)$', 'fontsize', fontsize, 'Interpreter','latex');
    set(gca, 'fontsize', fontsize);
    
    % Plot the input
    subplot(4,1,4);
    plot(tmat, u_mat, 'b', 'linewidth', 2);
    ylabel('u(t)');
    xlabel('time (sec)');
    set(gca, 'fontsize', fontsize);

end

function xdot = dynamics(t,x, u, P)
    % Extract the state and the estimate
    x_true = x(1:3);
    x_hat = x(4:6);
    
    % Calculate the state dynamics
    xdot_true = P.A * x_true + P.B*u;
    
    % Calculate the state estimate
    xdot_hat = P.A*x_hat + P.B*u + P.L*P.C*(x_true - x_hat);
    
    % Calculate the aggregate state dynamics
    xdot = [xdot_true; xdot_hat];
end

function P = controlParamters()
    % Define the desired values
    P.xd = zeros(3,1);
    
    % Define the system matrices
    P.A = [0 1 0; 0 -0.875 -20; 0 0 -50];
    P.B = [0; 0; 50];
    P.C = [1 0 0];
    
    % LQR definition
    gamma = 0.1;
    G = [1 0 0; 0 gamma 0];
    Q = G'*G;
    R = 0.01;
    P.K = lqr(P.A, P.B, Q, R);
    
    % LQG definition
    R = 1.0;
    Q = 10^8;
    Pkal = ss(P.A, P.B, P.C, 0);
    [~, P.L] = kalman(Pkal, inv(R), inv(Q));    
end
