function S04_L04_Linear_Helicopter()
%close all;

    %% Create the controller
    P = getSystemMatrices();
    %P = controlDesignPolePlacement(P);
    P = controlDesignLQR(P);
    
    % Set the control input function
    u = @(t, x) -P.K*x;
    
    %% Initialize the simulation variables
    % Time variables
    t0 = 0; % initial time
    dt = 0.1; % time step
    tf = 10.0; % final time
    t = t0:dt:tf;
    
    % set initial conditions:
    x0 = rand(8,1); % Random starting point
    xhat0 = rand(8,1); % Initial state estimate
    %xhat0 = x0;
    
    %% Simulate and plot the system using ode
    % Simulate the system           x       xhat          xhat
    [tvec xvec] = ode45(@(t,z) f(t,z(1:8), z(9:end), u(t,z(9:end)), P), t, [x0; xhat0]);
    tvec = tvec';
    xvec = xvec';  
    
    uvec = getControlVector(tvec, xvec(9:end, :), u);
    
    % Plot the resulting states
    plotResults(tvec, xvec(1:8,:), xvec(9:end, :), uvec);
end

function P = getSystemMatrices()
    %Rationalized Helicopter model
    %Code adapted from: http://folk.ntnu.no/skoge/book/2nd_edition/matlab_m/Sec13_2.m
    
    %% Create system matrices 
    % State matrix
    a01 = [          0                  0                  0   0.99857378005981;
                     0                  0   1.00000000000000  -0.00318221934140;
                     0                  0 -11.57049560546880  -2.54463768005371;
                     0                  0   0.43935656547546  -1.99818229675293;
                     0                  0  -2.04089546203613  -0.45899915695190;
    -32.10360717773440                  0  -0.50335502624512   2.29785919189453;
      0.10216116905212  32.05783081054690  -2.34721755981445  -0.50361156463623;
     -1.91097259521484   1.71382904052734  -0.00400543212891  -0.05741119384766];

    a02 = [0.05338427424431             0                  0                  0;
      0.05952465534210                  0                  0                  0;
     -0.06360262632370   0.10678052902222  -0.09491866827011   0.00710757449269;
                     0   0.01665188372135   0.01846204698086  -0.00118747074157;
     -0.73502779006958   0.01925575733185  -0.00459562242031   0.00212036073208;
                     0  -0.02121581137180  -0.02116791903973   0.01581159234047;
      0.83494758605957   0.02122657001019  -0.03787973523140   0.00035400385968;
                     0   0.01398963481188  -0.00090675335377  -0.29051351547241];

    P.A=[a01 a02];

    % Input matrix
    P.B=[              0                  0                  0                  0;
                      0                  0                  0                  0;
       0.12433505058289   0.08278584480286  -2.75247764587402  -0.01788876950741;
      -0.03635892271996   0.47509527206421   0.01429074257612                  0;
       0.30449151992798   0.01495801657438  -0.49651837348938  -0.20674192905426;
       0.28773546218872  -0.54450607299805  -0.01637935638428                  0;
      -0.01907348632812   0.01636743545532  -0.54453611373901   0.23484230041504;
      -4.82063293457031  -0.00038146972656                  0                 0];
  
    % Output matrix
    P.C = [1 0 0 0 0 0 0 0; ... % Pitch
           0 1 0 0 0 0 0 0; ... % Roll
           0 0 0 0 1 0 0 0; ... % Heading velocity (yaw rate)
           0 0 0 0 0 1 0 0]; ... % Heave velocity (forward velocity)
end

function P = controlDesignPolePlacement(P)
    %% Create feedback control matrix
    % Check controllability 
    rank_ctrb_mat = rank(ctrb(P.A, P.B))
    assert(rank_ctrb_mat == 8); % (must have rank of 8 or not completely controllable)
    
    % Create feedback control matrix
    lam = eig(P.A);
    %p = [lam(1:2); lam(5:end); -1.1; -1.2];
    p = [-1.1, -1.2, -1.3, -1.4, -1.5, -1.6, -1.7, -1.8];
    P.K = place(P.A, P.B, p);
    
    %% Create Luenberger observer matrix
    % Check observability
    rank_obsv_mat = rank(obsv(P.A, P.C))
    assert(rank_obsv_mat == 8); % (must have rank of 8 or not completely observable)
    
    % Create observability matrix
    p_o = [-10.1, -10.2, -10.3, -10.4, -10.5, -10.6, -10.7, -10.8].*10;
    P.L = place(P.A', P.C', p_o)'; % Note all of the transposes
end

function P = controlDesignLQR(P)
    %% Create feedback control matrix
    % Check controllability 
    rank_ctrb_mat = rank(ctrb(P.A, P.B))
    assert(rank_ctrb_mat == 8); % (must have rank of 8 or not completely controllable)
    
    % Create feedback control matrix
    Q = diag([1 1 1 1 1 1 1 1]);
    R = diag([10 10 10 10]).*10;
    P.K = lqr(P.A, P.B, Q, R);
    
    %% Create Luenberger observer matrix
    % Check observability
    rank_obsv_mat = rank(obsv(P.A, P.C))
    assert(rank_obsv_mat == 8); % (must have rank of 8 or not completely observable)
    
    % Create observability matrix
    Q = diag([100 1 1 1 1 1 1 1]);
    R = diag([1 1 1 1]);
    P.L = lqr(P.A', P.C', Q, R)'; % Note all of the transposes
end

function u_vec = getControlVector(tvec, xvec, u)
%getControlVector calculate the control vector over the specified time
%interval and state
%
% Inputs:
%   tvec: 1xm vector of time inputs
%   xvec: nxm matrix of states
%   u: function handle that takes time and state as inputs and outputs
%   the control input

    len = size(tvec, 2);
    u_vec = zeros(4, len);
    for k = 1:len
        u_vec(:,k) = u(tvec(k), xvec(:,k));
    end

end

function xdot = f(t, x, xhat, u, P)
    %f calculates the state dynamics using the current time, state, and
    %control input
    %
    % Inputs:
    %   t: current time
    %   x: current state
    %   xhat: current state estimate
    %   u: current control input
    %
    % Ouputs:
    %   xdot: time derivative of x(t)
    
    % LTI equation for state
    xtrue_dot = P.A*x + P.B*u;
    
    % LTI equation for state estimate
    y = P.C*x;
    xhat_dot = P.A*xhat + P.B*u + P.L*(y - P.C*xhat);
    
    xdot = [xtrue_dot; xhat_dot];
end

function plotResults(tvec, xvec, xhat_vec, uvec)

    %% Create figure and parameters
    % Figures: inputs and two figures for states
    fig_input = figure('units','normalized','outerposition',[0 0 1 1]);
    fig_state_2 = figure('units','normalized','outerposition',[0 0 1 1]);
    fig_state_1 = figure('units','normalized','outerposition',[0 0 1 1]);
    
    % Plot variables
    fontsize = 18;
    linewidth = 2;
    
    %% Plot desired value
    for k = 1:4
        set(0,'CurrentFigure', fig_state_1);
        subplot(4,1,k); hold on;
        plot([tvec(1) tvec(end)], [0 0], 'r:', 'linewidth', linewidth);
        
        set(0,'CurrentFigure', fig_state_2);
        subplot(4,1,k); hold on;
        plot([tvec(1) tvec(end)], [0 0], 'r:', 'linewidth', linewidth);
    end
    
    %% Plot the resulting state estimates
    linewidth = 3;
    color = 'k--';
    set(0,'CurrentFigure', fig_state_1);
    subplot(4,1,1); hold on;
    plot(tvec, xhat_vec(1,:), color, 'linewidth', linewidth);
    ylabel('$\theta$', 'fontsize', fontsize, 'Interpreter','latex');
    
    subplot(4,1,2); hold on;
    plot(tvec, xhat_vec(2,:), color, 'linewidth', linewidth);
    ylabel('$\phi$', 'fontsize', fontsize, 'Interpreter','latex');
    
    subplot(4,1,3); hold on;
    plot(tvec, xhat_vec(3,:), color, 'linewidth', linewidth);
    ylabel('$p$', 'fontsize', fontsize, 'Interpreter','latex');
    
    subplot(4,1,4); hold on;
    plot(tvec, xhat_vec(4,:), color, 'linewidth', linewidth);
    ylabel('$q$', 'fontsize', fontsize, 'Interpreter','latex');
    xlabel('Time (seconds)', 'fontsize', fontsize);
    
    set(0,'CurrentFigure', fig_state_2);
    subplot(4,1,1); hold on;
    plot(tvec, xhat_vec(5,:), color, 'linewidth', linewidth);
    ylabel('$\xi$', 'fontsize', fontsize, 'Interpreter','latex');
    
    subplot(4,1,2); hold on;
    plot(tvec, xhat_vec(6,:), color, 'linewidth', linewidth);
    ylabel('$v_x$', 'fontsize', fontsize, 'Interpreter','latex');
    
    subplot(4,1,3); hold on;
    plot(tvec, xhat_vec(7,:), color, 'linewidth', linewidth);
    ylabel('$v_y$', 'fontsize', fontsize, 'Interpreter','latex');
    
    subplot(4,1,4); hold on;
    plot(tvec, xhat_vec(8,:), color, 'linewidth', linewidth);
    ylabel('$v_z$', 'fontsize', fontsize, 'Interpreter','latex');
    xlabel('Time (seconds)', 'fontsize', fontsize);
    
    %% Plot the resulting states
    linewidth = 2;
    color = 'b';
    set(0,'CurrentFigure', fig_state_1);
    subplot(4,1,1); hold on;
    plot(tvec, xvec(1,:), color, 'linewidth', linewidth);
    ylabel('$\theta$', 'fontsize', fontsize, 'Interpreter','latex');
    
    subplot(4,1,2); hold on;
    plot(tvec, xvec(2,:), color, 'linewidth', linewidth);
    ylabel('$\phi$', 'fontsize', fontsize, 'Interpreter','latex');
    
    subplot(4,1,3); hold on;
    plot(tvec, xvec(3,:), color, 'linewidth', linewidth);
    ylabel('$p$', 'fontsize', fontsize, 'Interpreter','latex');
    
    subplot(4,1,4); hold on;
    plot(tvec, xvec(4,:), color, 'linewidth', linewidth);
    ylabel('$q$', 'fontsize', fontsize, 'Interpreter','latex');
    xlabel('Time (seconds)', 'fontsize', fontsize);
    
    set(0,'CurrentFigure', fig_state_2);
    subplot(4,1,1); hold on;
    plot(tvec, xvec(5,:), color, 'linewidth', linewidth);
    ylabel('$\xi$', 'fontsize', fontsize, 'Interpreter','latex');
    
    subplot(4,1,2); hold on;
    plot(tvec, xvec(6,:), color, 'linewidth', linewidth);
    ylabel('$v_x$', 'fontsize', fontsize, 'Interpreter','latex');
    
    subplot(4,1,3); hold on;
    plot(tvec, xvec(7,:), color, 'linewidth', linewidth);
    ylabel('$v_y$', 'fontsize', fontsize, 'Interpreter','latex');
    
    subplot(4,1,4); hold on;
    plot(tvec, xvec(8,:), color, 'linewidth', linewidth);
    ylabel('$v_z$', 'fontsize', fontsize, 'Interpreter','latex');
    xlabel('Time (seconds)', 'fontsize', fontsize);
    
    %% Plot inputs
    color = 'b';
    set(0,'CurrentFigure', fig_input);
    subplot(4,1,1); hold on;
    plot(tvec, uvec(1,:), color, 'linewidth', linewidth);
    ylabel('$u_1$', 'fontsize', fontsize, 'Interpreter','latex');
    
    subplot(4,1,2); hold on;
    plot(tvec, uvec(2,:), color, 'linewidth', linewidth);
    ylabel('$u_2$', 'fontsize', fontsize, 'Interpreter','latex');
    
    subplot(4,1,3); hold on;
    plot(tvec, uvec(3,:), color, 'linewidth', linewidth);
    ylabel('$u_3$', 'fontsize', fontsize, 'Interpreter','latex');
    
    subplot(4,1,4); hold on;
    plot(tvec, uvec(4,:), color, 'linewidth', linewidth);
    ylabel('$u_4$', 'fontsize', fontsize, 'Interpreter','latex');
    xlabel('Time (seconds)', 'fontsize', fontsize);
end

