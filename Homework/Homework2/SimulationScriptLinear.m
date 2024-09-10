function SimulationScriptLinear(A, B, dx0, xeq, du, ueq, title)
    %% Initialize the simulation variables
    t0 = 0; % initial time
    dt = 0.1; % time step
    tf = 10.0; % final time
        
    % Set the control input function
    % u = @(t, x) 0;  
    
    % Set the starting point
    % x0 = [0;0;0];

    %% Simulate and plot the system using ode
    % Simulate the system
    [tvec, dxvec] = matlabOde45(A, B, dx0, t0, dt, tf, du);
    duvec = getControlVector(tvec, dxvec, du);
    
    % Plot the resulting states
    figure;
    sgtitle(title);
    plotResults(tvec, dxvec, xeq, duvec, ueq, 'b');
    
    %% Simulate and plot the system using Euler (or other method)
    % % Simulate the system
    % [tvec, dxvec] = eulerIntegration(A, B, dx0, t0, dt, tf, du);
    % uvec = getControlVector(tvec, dxvec, ueq);
    % 
    % % Plot the results
    % plotResults(tvec, dxvec, xeq, duvec, ueq, 'r:');
    
    
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

function [tvec, xvec] = matlabOde45(A, B, x0, t0, dt, tf, u)
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
    odefunction = @(t,x) f(A,B,t,x,u(t,x));
    [tvec, xvec] = ode45(odefunction, t, x0);
    
    % Transpose the outputs to get in the correct form
    tvec = tvec';
    xvec = xvec';    
end

function [tvec, xvec] = eulerIntegration(A, B, x0, t0, dt, tf, u)
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
    % Write the euler simulation code here 
    % first state is initial state
    x = x0;
    for i = 2:len % from 2nd state to the end of time vector
        t = tvec(i-1); % get time from previous index (for passing into f)
        xdot = f(A, B, t, x, u(t, x)); % calculate first derivative at previous time, state, and input
        
        x = x + dt * xdot; % new state gets old state plus the scaled derivative at old state
        xvec(:,i) = x; % concatenate new state vector
    end
end

function zdot = f(A, B, t, z, u)
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
    % z2 = z(2);
    
    % Linear system matrices
    % A = [0,1;39.2,-156.8];
    % B = [1; 0];
    
    % LTI equation
    zdot = A*z + B*u;
end

function plotResults(tvec, dxvec, xeq, duvec, ueq, color)

    % Plot variables
    fontsize = 18;
    linewidth = 2;

    % Update with equilibrium offset
    xvec = dxvec + xeq;
    uvec = duvec + ueq;
    
    % Plot the resulting states
    subplot(3,1,1); hold on;
    plot(tvec, xvec(1,:), color, 'linewidth', linewidth);
    ylabel('theta(t)', 'fontsize', fontsize);
    
    subplot(3,1,2); hold on;
    plot(tvec, xvec(2,:), color, 'linewidth', linewidth);
    ylabel('thetadot(t)', 'fontsize', fontsize);
    
    subplot(3,1,3); hold on;
    plot(tvec, uvec, color, 'linewidth', linewidth);
    ylabel('u(t)', 'fontsize', fontsize);
    xlabel('time (s)', 'fontsize', fontsize);
end

