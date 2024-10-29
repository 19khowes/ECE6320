function S02_L03_PendulumEnergy(x0,xeq,ueq,K)
    % Define parameters
    P.g = 9.8; % Gravity constant
    P.l = 0.25; % Pendulum length
    P.m = 1/9.8; % Pendulum mass
    % P.b = 0.1; % Friction coefficient
    P.b = 1; % Friction coefficient

    % Simulate the state forward in time the state
    
    dt = 0.01;
    t = [0:dt:10];
    xd = xeq;
    u = @(t,x,P) ueq-K*(x-xeq);

    [tmat, xmat] = ode45(@(t,x)f(t,x,u, P), t, x0);
    tmat = tmat';
    xmat = xmat';
    umat = getControlVector(tmat, xmat, u);
    Tmat = min(1,max(-1,umat));
    
    % Calculate the energy
    len = length(tmat);
    E = zeros(1,len);
    for k = 1:len
        E(k) = calculateEnergy(xmat(:,k), P);
    end
    
    %% Plot the results
    fontsize = 12;
    
    % Plot the states
    figure;
    subplot(3,1,1);
    plot(tmat, xmat(1,:), 'b', 'linewidth', 3);
    ylabel('$\theta(t)$', 'fontsize', fontsize, 'Interpreter','latex');
    set(gca, 'fontsize', fontsize);
    
    subplot(3,1,2);
    plot(tmat, xmat(2,:), 'b', 'linewidth', 3);
    ylabel('$\dot{\theta}(t)$', 'fontsize', fontsize, 'Interpreter','latex');
    set(gca, 'fontsize', fontsize);

    % Plot the input
    subplot(3,1,3);
    plot(tmat, Tmat, 'b', 'linewidth', 3);
    ylabel('$T(t)$', 'fontsize', fontsize, 'Interpreter','latex');
    set(gca, 'fontsize', fontsize);
    
    % % Plot the energy
    % subplot(3,1,3);
    % plot(tmat, E, 'r', 'linewidth', 3);
    % ylabel('Total Energy');
    % xlabel('time (sec)');
    % set(gca, 'fontsize', fontsize);
    % 
    % %% Create a simulation
    % plotter = PendulumEnergyPlotter(P.l);
    % for k = 1:len
    %     plotter.plot(xmat(1,k), E(k), tmat(k));
    %     pause(dt/10);
    % end

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
        u_vec(:,k) = u(tvec(k), xvec(:,k), 0);
    end

end

function xdot = f(t,x, u_function, P)
    % Define parameters
    g = P.g;
    l = P.l;
    m = P.m;
    b = P.b;
    
    % Calculate the input
    u = u_function(t, x, P);
    
    % Threshold the control input
    T = max(-1, u);
    T = min(1, T);
    
    % Extract state
    theta = x(1);
    thetadot = x(2);
    
    % Calculate dynamics
    xdot = zeros(2,1);
    xdot(1) = thetadot;
    xdot(2) = g/l*sin(theta) - b/(m*l^2)*thetadot + 1/(m*l^2)*T;
end

function u = zeroControl(t, x, P)
    u = 0;
end

function e = calculateEnergy(x, P)
    % Define parameters
    g = P.g;
    l = P.l;
    m = P.m;
    b = P.b;
    
    % Extract state
    theta = x(1);
    thetadot = x(2);
    
    % Calculate position
    x_pos = l*sin(theta);
    y_pos = l*cos(theta);
    
    % Calculate velocity
    x_pos_dot = l*cos(theta)*thetadot;
    y_pos_dot = l*sin(theta)*thetadot;
    v = norm([x_pos_dot; y_pos_dot]);
    
    % Calculate the potential energy
    h = y_pos + l;
    PE = m*g*h;
    
    % Calculate the kinetic energy
    KE = 0.5*m*v^2;
    
    % Return total energy
    e = PE + KE;
end

