function S05_L02_RE_Numerical_Solution
close all;

    %% Initialize the sim
    % Create the time variables
    dt = 0.1; % Increase this by 0.1 from 0.1 to 0.5 
    t0 = 0;
    tf = 20;
    t_vec = t0:dt:tf;
    
    %% Collect the data
    % Get the data for the true P values
    [~, P_mat_true] = getTrueP(t_vec);
    
    % Get the data using euler integration
    [~, P_mat_euler] = getEulerP(t_vec, dt);
    
    % Get the data using matrix exponential
    [~, P_mat_exp] = getExpP(t_vec);
    
    
    %% Plot the P values
    % Create the figure
    fig_P = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    title('Solution values');
    
    % Plot true data
    linewidth = 3;
    fontsize = 12;
    color = 'b';    
    plotPmat(t_vec, P_mat_true, color, linewidth, fontsize);
    
    % Plot the euler data
    linewidth = 2;
    color = 'r';
    plotPmat(t_vec, P_mat_euler, color, linewidth, fontsize);
    
    % Plot the exponential data
    linewidth = 2;
    color = 'g:';
    plotPmat(t_vec, P_mat_exp, color, linewidth, fontsize);
    
    %% Simulate the state forward in time
    % Set initial conditions
    x0 = [5; 5];
    
    % Simulate the true response
    [~, x_mat_true] = simulateODE(t_vec, x0);
    
    % Simulate using Euler data
    [~, x_mat_euler] = simulateEuler(t_vec, dt, P_mat_euler, x0);
    
    % Simulate using exponential data
    [~, x_mat_exp] = simulateEuler(t_vec, dt, P_mat_exp, x0);
    
    %% Plot the states
    figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    title('State Simulation');
    
    % Plot true data
    linewidth = 3;
    fontsize = 12;
    color = 'b';
    plotStates(t_vec, x_mat_true, color, linewidth, fontsize);
    
    % Plot the Euler data
    linewidth = 2;
    color = 'r';
    plotStates(t_vec, x_mat_euler, color, linewidth, fontsize);
    
    % Plot the exponential data
    linewidth = 2;
    color = 'g:';
    plotStates(t_vec, x_mat_exp, color, linewidth, fontsize);
    
    % Set the P plot in front
    figure(fig_P);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotPmat(t_vec, Pmat, color, linewidth, fontsize)
    [p11, p12, p21, p22] = extractRows(Pmat);
    
    % Plot first row, first column
    subplot(4,1,1); hold on;
    plot(t_vec, p11, color, 'linewidth', linewidth);
    ylabel('p_{11}', 'fontsize', fontsize);
    
    % Plot first row, second column
    subplot(4,1,2); hold on;
    plot(t_vec, p12, color, 'linewidth', linewidth);
    ylabel('p_{12}', 'fontsize', fontsize);
    
    % Plot second row, first column
    subplot(4,1,3); hold on;
    plot(t_vec, p21, color, 'linewidth', linewidth);
    ylabel('p_{21}', 'fontsize', fontsize);
    
    % Plot second row, second column
    subplot(4,1,4); hold on;
    plot(t_vec, p22, color, 'linewidth', linewidth);
    ylabel('p_{22}', 'fontsize', fontsize);
    xlabel('time');
end

function [p11, p12, p21, p22] = extractRows(P_mat)
    p11 = squeeze(P_mat(1,1,:));
    p12 = squeeze(P_mat(1,2,:));
    p21 = squeeze(P_mat(2,1,:));
    p22 = squeeze(P_mat(2,2,:));
end

function plotStates(t_vec, x_mat, color, linewidth, fontsize)
    % Plot first row, first column
    subplot(2,1,1); hold on;
    plot(t_vec, x_mat(1,:), color, 'linewidth', linewidth);
    ylabel('x_{1}', 'fontsize', fontsize);
    
    % Plot first row, second column
    subplot(2,1,2); hold on;
    plot(t_vec, x_mat(2,:), color, 'linewidth', linewidth);
    ylabel('x_{2}', 'fontsize', fontsize);
    xlabel('time');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% True solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t_vec, P_mat] = getTrueP(t_vec)
%getTrueP Loops through the time input and calculates the P matrix

    % Initialize the P_mat for all time indices
    len = length(t_vec);
    P_mat = zeros(2,2,len);
    
    % Loop through and get the matrix
    T = t_vec(end);
    for k = 1:len
        P_mat(:,:,k) = trueSolution(t_vec(k), T);
    end    
end

function P = trueSolution(t, T)
%trueSolution: Provides the true solution to the Riccati equation at time t
%given the final time T
    P = 1/(1 + 1/3*(T-t)^3)*[1 T-t; T-t (T-t)^2];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Euler solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t_vec, P_mat] = getEulerP(t_vec, dt)
    % Initialize P_mat
    len = length(t_vec);
    P_mat = zeros(2,2,len);
    P_mat(:,:,len) = [1 0; 0 0]; % Final P value = S
    
    % Initialize state matrices
    A = [0 1; 0 0];
    B = [0; 1];
    
    % Initialize costs
    Q = zeros(2);
    R = 1;
    R_inv = 1/R;
    
    % move from end to beginning
    P = squeeze(P_mat(:,:,end));
    for k = len:-1:2
        % Calculate the time derivative
        P_dot = -A'*P - P*A - Q + P*B*R_inv*B'*P;
        
        % Update the state moving backward in time
        P = P - dt * P_dot;
        
        % Store the state
        P_mat(:,:,k-1) = P;
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matrix exponential solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t_vec, P_mat] = getExpP(t_vec)
    % Initialize P_mat
    len = length(t_vec);
    P_mat = zeros(2,2,len);
    
    % Initialize state matrices
    A = [0 1; 0 0];
    B = [0; 1];
    
    % Initialize costs
    Q = zeros(2);
    R = 1;
    S = [1 0; 0 0]; % Final P value = S        
    
    % Set initial conditions for X and Y
    X0 = eye(2);
    Y0 = S;
    
    % Create aggregate matrix
    M = [A -B*inv(R)*B'; -Q -A'];
    
    % Loop through and calculate P
    T = t_vec(end);
    for k = 1:len
        % Calculate X and Y
        t = t_vec(k);
        agg_mat = expm(M*(t-T)) * [X0; Y0];
        X = agg_mat(1:2, :);
        Y = agg_mat(3:4, :);
        
        % Calculate P(t)
        P_mat(:,:, k) = Y/X;
    end
end

function P = calculateP(t, T)
    % Initialize state matrices
    A = [0 1; 0 0];
    B = [0; 1];
    
    % Initialize costs
    Q = zeros(2);
    R = 1;
    S = [1 0; 0 0]; % Final P value = S        
    
    % Set initial conditions for X and Y
    X0 = eye(2);
    Y0 = S;
    
    % Create aggregate matrix
    M = [A -B*R*B'; -Q -A'];
    
    agg_mat = expm(M*(t-T)) * [X0; Y0];
    X = agg_mat(1:2, :);
    Y = agg_mat(3:4, :);
    
    % Calculate P
    P = Y/X;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t_vec, x_mat] = simulateEuler(t_vec, dt, P_mat, x0)
    % Initialize x_mat
    len = length(t_vec);
    x_mat = zeros(2,len);
    
    % Initialize state matrices
    A = [0 1; 0 0];
    B = [0; 1];
    
    % Initialize costs
    Q = zeros(2);
    R = 1;
    R_inv = 1/R;
    
    % Set initial state
    x = x0;
    x_mat(:,1) = x0;
    
    % Loop forward in time using euler integration
    for k = 1:len-1
        % Calculate control
        P = squeeze(P_mat(:,:,k));
        u = -R_inv * B'*P*x;
        
        % Calculate time derivative
        xdot = A*x + B*u;
        
        % Update and store state
        x = x + dt*xdot;
        x_mat(:,k+1) = x;
    end
end

function [t_vec, x_mat] = simulateODE(t_vec, x0)
    T = t_vec(end);
    [~, x_mat] = ode45(@(t,x)dynamics(t, x, T), t_vec, x0);
    x_mat = x_mat';
end

function xdot = dynamics(t, x, T)
    % Initialize state matrices
    A = [0 1; 0 0];
    B = [0; 1];
    
    % Initialize costs
    Q = zeros(2);
    R = 1;
    R_inv = 1/R;
    
    % Calculate P
    P = calculateP(t, T);
    
    % Calculate control input
    u = -R_inv*B'*P*x;
    
    % Calculate dynamics
    xdot = A*x + B*u;
end