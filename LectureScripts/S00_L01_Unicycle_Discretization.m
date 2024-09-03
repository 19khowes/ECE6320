function S00_L01_Unicycle_Discretization
close all;

    %% Initialize the simulation variables
    % Timing variables
    t0 = 0;     % Initial time
    dt = 1;   % Step size
    tf = 10.0;  % Final time
    
    % Initial state
    x0 = [0;0;0];
    
    % Set the control input function
    u = @(t, x) [1; 0.25];
    %u = @(t, x) [sin(t); cos(t)];
    %u = @(t, x) [t^2*sin(10*t); 1+3*cos(4*t)];
    
    %% Simulate the system using ode45 - Continuous
    % Simulate the system
    [tvec, xvec_true] = ode45(@(t,x) unicycleDynamics(t, x, u(t,x)), t0:dt:tf, x0);
    xvec_true = xvec_true'; % Make each state a column vector
    uvec = getControlVector(tvec, xvec_true, u);
    
    % Plot the resulting states
    figure('units','normalized','outerposition',[0 0 1 1]);
    plotResults(tvec, xvec_true, uvec, 'b');
    
    %% Simulate the system using ode45 - discrete
    % Simulate the system
    [tvec, xvec_true] = matlabOde45Discrete(x0, t0, dt, tf, u, @unicycleDynamics);
    uvec = getControlVector(tvec, xvec_true, u);
    
    % Plot the resulting states
    %figure('units','normalized','outerposition',[0 0 1 1]);
    plotResults(tvec, xvec_true, uvec, 'b');
    
    %% Simulate the system using an Euler discretization
    % Simulate the system
    [tvec, xvec_euler] = eulerIntegration(x0, t0, dt, tf, u, @unicycleDynamics);
    uvec = getControlVector(tvec, xvec_euler, u);
    
    % Plot the resulting states
    plotResults(tvec, xvec_euler, uvec, 'ro');
    
    % Plot the resulting error
    plotError(tvec, xvec_true, xvec_euler, 'r');
    
    %% Simulate the system using an Rk4 discretization
    % Simulate the system
    [tvec, xvec_rk4] = rk4Integration(x0, t0, dt, tf, u, @unicycleDynamics);
    uvec = getControlVector(tvec, xvec_rk4, u);
    
    % Plot the resulting states
    plotResults(tvec, xvec_rk4, uvec, 'go');
    
    % Plot the resulting error
    plotError(tvec, xvec_true, xvec_rk4, 'g');
    
end

function xdot = unicycleDynamics(~, x, u)
% unicycleDynamics Calculates the dynamics of the unicycle robot given the
% inputs
%
% Inputs:
%   t: Current time (not used)
%   x: Unicycle configuration ( (x,y) position and orientation)
%   u: Input to the unicycle
%       u(1): Translational velocity input
%       u(2): Rotational velocity input
%
% Outputs:
%   xdot: dynamics at time t

    % Extract orientation
    th = x(3);
    
    % Extract inputs
    v = u(1);
    w = u(2);
    
    % Calculate the dynamics
    xdot = [ v * cos(th); ...
             v * sin(th); ...
             w];
end

function [tvec, xvec] = matlabOde45Discrete(x0, t0, dt, tf, u, f)
    %MatlabOde45 uses ODE 45 to simulate the state starting at x0 from time
    % t0 to tf assuming discrete steps in the input (i.e., the input is
    % only updated every dt seconds)
    %
    % Inputs:
    %   x0: nx1 initial state
    %   t0: scalar - initial time
    %   dt: scalar - time increment
    %   tf: scalar - final time
    %   u: function handle that takes time and state as inputs and outputs
    %   the control input
    %   f: function handle that take in time, state, and input
    %
    % Outputs:
    %   tvec: 1xm vector of times associated with the states
    %   xvec: nxm matrix of states where each column is a state at the
    %   associated time in tvec
    
    % Initialize values
    tvec = t0:dt:tf;
    len = length(tvec);
    xvec = zeros(size(x0,1), len);
    xvec(:,1) = x0;
    
    % Simulate forward in time
    x = x0;
    for k = 2:len
        % Calculate the control input
        t = tvec(k-1);  % Time at step k-1
        u_km1 = u(t,x); % Input at step k-1
        
        % Simulate forward in time
        [~, xmat] = ode45(@(t,x) f(t, x, u_km1), [t t+dt], x);
        x = xmat(end,:)'; % Extract the final state from xmat
        xvec(:,k) = x;
    end
end

function [tvec, xvec] = eulerIntegration(x0, t0, dt, tf, u, f)
    %eulerIntegration uses eulerIntegration to simulate the state starting at x0 from time
    % t0 to tf
    %
    % Inputs:
    %   x0: nx1 initial state
    %   t0: scalar - initial time
    %   dt: scalar - time increment
    %   tf: scalar - final time
    %   u: function handle that takes time and state as inputs and outputs
    %   the control input
    %   f: function handle that take in time, state, and input
    %
    % Outputs:
    %   tvec: 1xm vector of times associated with the states
    %   xvec: nxm matrix of states where each column is a state at the
    %   associated time in tvec
    
    % Initialize values
    tvec = t0:dt:tf;
    len = length(tvec);
    xvec = zeros(size(x0,1), len);
    xvec(:,1) = x0;
    
    % Simulate forward in time
    x = x0;
    for k = 2:len
        % Calculate state dynamics
        t = tvec(k-1);
        xdot = f(t, x, u(t,x));
        
        % Simulate forward in time
        x = x + dt * xdot;
        xvec(:,k) = x;
    end
end

function [tvec, xvec] = rk4Integration(x0, t0, dt, tf, u, f)
    %rk4Integration uses a four step runge-kutta discretization of the system
    %
    % Inputs:
    %   x0: nx1 initial state
    %   t0: scalar - initial time
    %   dt: scalar - time increment
    %   tf: scalar - final time
    %   u: function handle that takes time and state as inputs and outputs
    %   the control input
    %   f: function handle that take in time, state, and input
    %
    % Outputs:
    %   tvec: 1xm vector of times associated with the states
    %   xvec: nxm matrix of states where each column is a state at the
    %   associated time in tvec
    
    % Initialize values
    tvec = t0:dt:tf;
    len = length(tvec);
    xvec = zeros(size(x0,1), len);
    xvec(:,1) = x0;
    
    % Simulate forward in time
    x = x0;
    for k = 2:len
        % Calculate state dynamics
        t = tvec(k-1);
        u_km1 = u(t,x); % Input at step k-1
        k1 = f(t, x, u_km1);
        k2 = f(t+dt/2, x+(dt/2)*k1, u_km1);
        k3 = f(t+dt/2, x+(dt/2)*k2, u_km1);
        k4 = f(t+dt, x+dt*k3, u_km1);
        
        % Simulate forward in time
        x = x + (dt/6)*k1 + (dt/3)*k2 + (dt/3)*k3 + (dt/6)*k4;
        xvec(:,k) = x;
    end
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
    u_vec = zeros(2, len);
    for k = 1:len
        u_vec(:,k) = u(tvec(k), xvec(:,k));
    end

end

function plotResults(tvec, xvec, uvec, color)

    % Plot variables
    fontsize = 18;
    linewidth = 2;
    
    % Plot the resulting states
    subplot(3,2,1); hold on;
    plot(tvec, xvec(1,:), color, 'linewidth', linewidth);
    ylabel('q_1(t)', 'fontsize', fontsize);
    
    subplot(3,2,3); hold on;
    plot(tvec, xvec(2,:), color, 'linewidth', linewidth);
    ylabel('q_2(t)', 'fontsize', fontsize);
    
    subplot(3,2,5); hold on;
    plot(tvec, xvec(3,:), color, 'linewidth', linewidth);
    ylabel('q_3(t)', 'fontsize', fontsize);
    
%     subplot(5,1,4); hold on;
%     plot(tvec, uvec(1,:), color, 'linewidth', linewidth);
%     ylabel('v(t)', 'fontsize', fontsize);
%     xlabel('time (s)', 'fontsize', fontsize);
%     
%     subplot(5, 1, 5); hold on;
%     plot(tvec, uvec(2,:), color, 'linewidth', linewidth);
%     ylabel('\omega(t)', 'fontsize', fontsize);
%     xlabel('time (s)', 'fontsize', fontsize);
end

function plotError(tvec, xvec_truth, xvec, color)
    % Plot variables
    fontsize = 18;
    linewidth = 2;
    
    % Calculate the normalized difference
    x_err = abs(xvec_truth - xvec);
    
    % Plot the resulting states
    subplot(3,2,2); hold on;
    plot(tvec, x_err(1,:), color, 'linewidth', linewidth);
    ylabel('q_1(t)', 'fontsize', fontsize);
    
    subplot(3,2,4); hold on;
    plot(tvec, x_err(2,:), color, 'linewidth', linewidth);
    ylabel('q_2(t)', 'fontsize', fontsize);
    
    subplot(3,2,6); hold on;
    plot(tvec, x_err(3,:), color, 'linewidth', linewidth);
    ylabel('q_3(t)', 'fontsize', fontsize);

end