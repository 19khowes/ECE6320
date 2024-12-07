function thermostat_control_simulation()
    close all;
    clear all;

    %% Create system
    [A, B, C, E] = get_system();
    n = height(A);
    m = width(B);
        
    %% Control Design
    syms x1_d x2_d x4_d x5_d uff1 uff2;
    x_d = [x1_d; x2_d; 20; x4_d; x5_d];
    xddot = [0;0;0;0;0];
    % uff = [uff1; uff2];
    uff = [uff1; 0]; % uff2 is free, choose 0
    eqn = B*uff + A*x_d - xddot == 0; % solve for uff and x_d (excpet x_d3)
    sol = solve(eqn);
    uff = double(subs(uff, uff1, sol.uff1));
    x_d = double(subs(x_d, x_d, [sol.x1_d; sol.x2_d; 20; sol.x4_d; sol.x5_d]));
    % Sanity check
    check = B*uff + A*x_d - xddot;
    % check controllability
    Gamma = ctrb(A,B);
    rGamma = rank(Gamma); % fully controllable
    % LQR
        % form aggregate state
    C1 = [0 0 1 0 0];
    C2 = [1];
    Abar = [A zeros([n 2]); C1 0 0; 0 0 0 0 0 C2 0];
    Bbar = [B; zeros([2 m])];
        % pick and verify Q/R
    Q = diag([0 0 1 0 0 (1/(20^2)) (1/(20^2))]);
    R = diag([(1/(0.5^2)) (1/(0.5^2))]);
    Cc = sqrt(Q);
    checkQ = Cc'*Cc;
    Omegac = obsv(Abar, Cc);
    rOmegac = rank(Omegac);
    % Create K using lqr()
    K = lqr(Abar,Bbar,Q,R);
    % check eigenvalues of closed loop system
    eigclosed = eig(Abar-Bbar*K);


    %% Create the observer (Create the observer)
    T = [eye(5); zeros([2 5])];
    Kbar = K*T;
    % check observability
    Omega = obsv(A,C);
    rOmega = rank(Omega); % completely observable
    QL = diag(10*ones([1 n]));
    RL = 1;
    L = lqr(A',C',QL,RL)';
    % check observer convergence
    eigALC = eig(A-L*C);

    
    
    %% Store the values (Feel free to add any additional values here to pass into the control or dynamics functions)
    P.A = A;
    P.B = B;
    P.C = C;
    P.E = E;
    P.n_states = 5;
    P.n_ctrl = 2;
    P.ctrl_max = 100;   

    P.uff = uff;
    P.x_d = x_d;
    P.K = K;
    P.L = L;
    P.d_nom = 10;
    

    %% Simulate the system (Don't change anything in this section except x0_ctrl)
    % Create the initial state
    x0_sys = [20.5; 19.; 19.; 11.; 21];
    x0_obs = [15; 15; 15; 15; 15];
    x0_ctrl = [0; 0]; % Any additional states added for control -> initialize each to zero
    x0 = [x0_sys; x0_obs; x0_ctrl];

    % Simulate the system throughout the entire day (86400 seconds)
    t_int = [0 86400];
    [tvec, x_mat] = ode45(@(t,x)dynamics(t, x, P), t_int, x0);
    x_mat = x_mat';

    % Get the control and disturbance over time
    u_mat = get_all_control(x_mat, P);
    d_vec = get_all_disturbances(tvec);
    

    %% Plot the results
    % Plot the states over time
    figure;
    for k = 1:5
        subplot(5,1,k);
        plot([tvec(1) tvec(end)]/3600, [x_d(k) x_d(k)], 'r:')        
        hold on;
        plot(tvec/3600, x_mat(k+P.n_states, :), 'g')
        plot(tvec/3600, x_mat(k,:), 'b');
        ylabel(["x_", num2str(k)])
    end
    xlabel("Time (hr)")
    sgtitle("States vs Time")  
    legend('Desired', 'Observed States', 'States');

    % Plot other things vs time
    % Error vs time
    figure;
    for k = 1:5
        subplot(5,1,k);
        plot(tvec/3600, x_mat(k,:)-x_d(k), 'b');
        ylabel(["x_", num2str(k), " Error"])
    end
    xlabel("Time (hr)")
    sgtitle("State Error vs Time") 

    % Control vs time
    figure;
    for k = 1:2
        subplot(2,1,k);
        plot(tvec/3600, u_mat(k,:), 'b');
        ylabel(["u_", num2str(k)])
    end
    xlabel("Time (hr)")
    sgtitle("Control vs Time") 

    % Disturbance vs time
    figure;
    plot(tvec/3600, d_vec, 'b');
    ylabel("Disturbance")
    xlabel("Time (hr)")
    sgtitle("Disturbance vs Time") 

end

function [A, B, C, E] = get_system()
    data = load("sys_1.mat");
    A = data.A;
    B = data.B;
    C = data.C;
    E = data.E;
end

function d = disturbance(t)
    % Disturbance as a function of time
    if t < 3600 || t > 82800
        d = 9.0;
    elseif t < 43200
        delta = (t - 3600)/(43200-3600);
        d = 15*(delta) + 9*(1-delta);
    elseif t < 82800
        delta = (t - 43200)/(82800-43200);
        d = 9*(delta) + 15*(1-delta);
    end
    
end

function xdot = dynamics(t, x, P)
    % Defines the dynamics of the system
    % Args:
    %   t: Current time
    %   x: Current augmented state (i.e., actual state, observer state, and
    %   states needed for control)
    %   P: Parameters that are required for simulation
    %
    % Returns:
    %   xdot: The aggregate dynamics

    % Extract the states (don't touch the next three lines)
    x_sys = x(1:P.n_states, :);
    x_obs = x(P.n_states+1:2*P.n_states, :);
    x_ctrl = x(2*P.n_states+1:end, :);

    % Calculate the control (don't touch this line)
    u = control(x_obs, x_ctrl, P);

    % System dynamics (don't touch this line)
    x_sys_dot = P.A*x_sys + P.B*u + P.E*disturbance(t);

    % Observer dynamics (this line should use P.d_nom instead of disturbance(t))
    y = P.C*x_sys; % sensor measurement (output)
    x_obs_dot = P.A*x_obs + P.B*u + P.L*(y - P.C*x_obs) + P.E*P.d_nom; % TODO: should this have a d value in it?

    % Control dynamics (definitely fix this line)
    x_ctrl_dot = [x_sys(3) - P.x_d(3); x_ctrl(1)];
    % x_ctrl_dot = [0;0];

    xdot = [x_sys_dot; x_obs_dot; x_ctrl_dot];
end

function u = control(x_obs, x_ctrl, P)
    % Defines the control 
    % 
    % Args:
    %   x_obs: The observer state
    %   x_ctrl: The controller state
    %   P: Parameters that are required for control
    %
    % Returns
    %   u: The resulting control

    % Create the feedback control (you'll want to change this)
    % u = zeros(P.n_ctrl,1);
    u = P.uff - P.K*([x_obs-P.x_d; x_ctrl]);

    % Bound the control (leave the following two lines alone)
    u = max(u, -P.ctrl_max);
    u = min(u, P.ctrl_max);
end

function u_mat = get_all_control(x_mat, P)
    % Loop through and calculate the control given the states
    %
    % Args:
    %   x_mat: matrix of all states over time
    %   P: Necessary parameters for computing control

    % Extract the observer and control states
    x_mat_obs = x_mat(P.n_states+1:2*P.n_states, :);
    x_mat_ctrl = x_mat(2*P.n_states+1:end, :);
    
    % Calculate the control to be returned
    n_times = size(x_mat, 2);
    u_mat = zeros(P.n_ctrl, n_times);
    for k = 1:n_times
        u_mat(:,k) = control(x_mat_obs(:,k), x_mat_ctrl(:,k), P);
    end    
end

function d_vec = get_all_disturbances(tvec)
    % Loop through and calculates the disturbance given the time values
    n_times = length(tvec);
    d_vec = zeros(1,n_times);
    for k = 1:n_times
        d_vec(k) = disturbance(tvec(k));
    end
end
