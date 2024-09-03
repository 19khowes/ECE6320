function S03_L04_LQR_PendulumCart_Solution()
    %close all;
    
    %% Initialize the simulation variables
    t0 = 0; % initial time
    dt = 0.01; % time step
    tf = 10.0; % final time
    x0 = [0;0;pi-.25; 0];
    
    % Calculate the feedback matrix
    %K_part2 = calculateFeedbackMatrixPart2();
    K_part2 = calculateFeedbackMatrixPart2_LQR();
    
    % Set the control input function
    u = @(t, x) input_part_2(t, x, K_part2);   
        
    %% Simulate and plot the system using ode
    % Simulate the system
    [tvec, xvec] = matlabOde45(x0, t0, dt, tf, u);
    uvec = getControlVector(tvec, xvec, u);
    
    % Plot the resulting states
    figure;
    plotResults(tvec, xvec, uvec, 'b');
end

function K = calculateFeedbackMatrixPart2()
    % create the feedback control matrix for part 2 (linearized about theta = pi)
    A = [0 1 0 0; 0 -10/55 147/55 0; 0 0 0 1; 0 -5/11 343/11 0]; 
    B = [0; 100/55; 0; 50/11];
    %eig_A = eig(A); % Calculate eigenvalues, only change the positive values
    P = [-0.14283 -0.6041 -1 -10];
    K = acker(A, B, P) % Create the feedback control matrix
    
    eig_A_BK = eig(A-B*K)
end

function K = calculateFeedbackMatrixPart2_LQR()
    % create the feedback control matrix for part 2 (linearized about theta = pi)
    A = [0 1 0 0; 0 -10/55 147/55 0; 0 0 0 1; 0 -5/11 343/11 0]; 
    B = [0; 100/55; 0; 50/11];
    
    % Create cost matrices
    Q = diag([10 .1 .1 100]);
    R = 1;
    
    K = lqr(A, B, Q, R) % Create the feedback control matrix
    
    eig_A_BK = eig(A-B*K)
end

function u = input_part_2(t, z, K)
%input_part_2 calculates the input for homework 3 part 2 where the system
%is linearized about theta = pi
%
% Inputs:
%   t: current time (not used)
%   z: current state [x; xdot; theta; theta_dot]
%   K: feedback matrix

    % Calculate dz
    z_eq = [0; 0; pi; 0];
    dz = z - z_eq;
    
    % Calculate the control
    u = -K*dz;
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
    u_vec = zeros(1, len);
    for k = 1:len
        u_vec(:,k) = u(tvec(k), xvec(:,k));
    end

end

function [tvec, xvec] = matlabOde45(x0, t0, dt, tf, u)
    %MatlabOde45 uses ODE 45 to simulate the state starting at x0 from time
    % t0 to tf
    %
    % Inputs:
    %   x0: nx1 initial state
    %   t0: scalar - initial time
    %   dt: scalar - time increment
    %   tf: scalar - final time
    %   u: function handle that takes time and state as inputs and outputs
    %   the control input
    %
    % Outputs:
    %   tvec: 1xm vector of times associated with the states
    %   xvec: nxm matrix of states where each column is a state at the
    %   associated time in tvec
    
    % Initialize the time
    t = t0:dt:tf;
    
    % Simulate the output
    [tvec xvec] = ode45(@(t,x) f(t,x,u(t,x)), t, x0);
    
    % Transpose the outputs to get in the correct form
    tvec = tvec';
    xvec = xvec';    
end

function zdot = f(t, z, u)
    %f calculates the state dynamics using the current time, state, and
    %control input
    %
    % Inputs:
    %   t: current time
    %   x: current state
    %   u: current control input
    %
    % Ouputs:
    %   xdot: time derivative of x(t)
    
    % Extract states
    xdot = z(2); % xdot - the time derivative of the position
    theta = z(3); % theta - angle of the pendulum from verticle
    thetadot = z(4); % thetadot - the time derivative of the angle
    
    % Calculate the second derivatives
    xddot = Xddot(thetadot, theta, u, xdot);
    thetaddot = Thetaddot(thetadot, theta, u, xdot);
    
    % Output the dynamics
    zdot = [xdot; xddot; thetadot; thetaddot];
end

function plotResults(tvec, xvec, uvec, color)

    % Plot variables
    fontsize = 18;
    linewidth = 2;
    
    % Plot the resulting states
    subplot(5,1,1); hold on;
    plot(tvec, xvec(1,:), color, 'linewidth', linewidth);
    ylabel('$x(t)$', 'fontsize', fontsize, 'Interpreter','latex');
    
    subplot(5,1,2); hold on;
    plot(tvec, xvec(2,:), color, 'linewidth', linewidth);
    ylabel('$\dot{x}(t)$', 'fontsize', fontsize, 'Interpreter','latex');
    
    subplot(5,1,3); hold on;
    plot(tvec, xvec(3,:), color, 'linewidth', linewidth);
    ylabel('$\theta(t)$', 'fontsize', fontsize, 'Interpreter','latex');
    
    subplot(5,1,4); hold on;
    plot(tvec, xvec(4,:), color, 'linewidth', linewidth);
    ylabel('$\dot{\theta}(t)$', 'fontsize', fontsize, 'Interpreter','latex');
    
    subplot(5,1,5); hold on;
    plot(tvec, uvec, color, 'linewidth', linewidth);
    ylabel('$u(t)$', 'fontsize', fontsize, 'Interpreter','latex');
    xlabel('Time (seconds)', 'fontsize', fontsize);
    
    % Loop through and plot the result
    figure;
    cart = CartPendulumPlotter();
    len = length(tvec);
    for k = 1:len
        t = tvec(k);
        cart.plot(xvec(1,k), xvec(3,k));
        pause(0.01);
    end
end

function xddot = Xddot(thetadot,theta,u,xdot)
%XDDOT
%    XDDOT = XDDOT(THETADOT,THETA,U,XDOT)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    11-Sep-2019 21:13:46

t2 = cos(theta);
t3 = sin(theta);
xddot = ((u.*1.0e2-xdot.*1.0e1+t2.*t3.*1.47e2+t3.*thetadot.^2.*6.0).*(-1.0./5.0))./(t2.^2.*3.0-1.4e1);

end

function thetaddot = Thetaddot(thetadot,theta,u,xdot)
%THETADDOT
%    THETADDOT = THETADDOT(THETADOT,THETA,U,XDOT)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    11-Sep-2019 21:13:45

t2 = cos(theta);
t3 = sin(theta);
thetaddot = (t3.*3.43e2+t2.*u.*5.0e1-t2.*xdot.*5.0+t2.*t3.*thetadot.^2.*3.0)./(t2.^2.*3.0-1.4e1);

end

