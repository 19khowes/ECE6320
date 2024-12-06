function thermostat_control_simulation()
    close all;
    clear all;

    %% Create system
    [A, B, C, E] = get_system();
        
    %% Control Design
    x_d = [0.; 0.; 20.; 0.; 0.]; % This is definitely not correct


    %% Create the observer (Create the observer)
    
    
    %% Store the values (Feel free to add any additional values here to pass into the control or dynamics functions)
    P.A = A;
    P.B = B;
    P.C = C;
    P.E = E;
    P.n_states = 5;
    P.n_ctrl = 2;
    P.ctrl_max = 100;    
    

    %% Simulate the system (Don't change anything in this section except x0_ctrl)
    % Create the initial state
    x0_sys = [20.5; 19.; 19.; 11.; 21];
    x0_obs = [15; 15; 15; 15; 15];
    x0_ctrl = []; % Any additional states added for control -> initialize each to zero
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

    % Plot other things vs time

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
    x_obs_dot = zeros(size(x_obs));

    % Control dynamics (definitely fix this line)
    x_ctrl_dot = zeros(size(x_ctrl));

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
    u = zeros(P.n_ctrl,1);

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
