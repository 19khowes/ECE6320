function PlotCart(z0,zeq,ueq,K)
    % Define parameters

    % Simulate the state forward in time the state
    
    dt = 0.01;
    t = [0:dt:10];
    u = @(t,z) -K*(z-zeq(t,z))+ueq(t,z);

    [tmat, xmat] = ode45(@(t,x)f(t,x,u), t, z0);
    tmat = tmat';
    xmat = xmat';
    umat = getControlVector(tmat, xmat, u);
    
    %% Plot the results
    fontsize = 12;
    
    % Plot the states
    figure;
    subplot(5,1,1);
    plot(tmat, xmat(1,:), 'b', 'linewidth', 3);
    ylabel('$x(t)$', 'fontsize', fontsize, 'Interpreter','latex');
    set(gca, 'fontsize', fontsize);
    
    subplot(5,1,2);
    plot(tmat, xmat(2,:), 'b', 'linewidth', 3);
    ylabel('$\dot{x}(t)$', 'fontsize', fontsize, 'Interpreter','latex');
    set(gca, 'fontsize', fontsize);

    subplot(5,1,3);
    plot(tmat, xmat(3,:), 'b', 'linewidth', 3);
    ylabel('$\theta(t)$', 'fontsize', fontsize, 'Interpreter','latex');
    set(gca, 'fontsize', fontsize);

    subplot(5,1,4);
    plot(tmat, xmat(4,:), 'b', 'linewidth', 3);
    ylabel('$\dot{\theta}(t)$', 'fontsize', fontsize, 'Interpreter','latex');
    set(gca, 'fontsize', fontsize);

    % Plot the input
    subplot(5,1,5);
    plot(tmat, umat, 'b', 'linewidth', 3);
    ylabel('$u(t)$', 'fontsize', fontsize, 'Interpreter','latex');
    set(gca, 'fontsize', fontsize);
    

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

function zdot = f(t,z, u_function)
    % Define parameters
    
    % Calculate the input
    u = u_function(t, z);
    
    % Extract state
    x = z(1);
    xdot = z(2);
    theta = z(3);
    thetadot = z(4);
    
    % Calculate dynamics
    zdot = zeros(4,1);
    zdot(1) = xdot;
    zdot(2) = -(6*sin(theta)*thetadot^2 + 100*u - 10*xdot + 147*cos(theta)*sin(theta))/(5*(3*cos(theta)^2 - 14));
    zdot(3) = thetadot;
    zdot(4) = (3*cos(theta)*sin(theta)*thetadot^2 + 343*sin(theta) + 50*u*cos(theta) - 5*xdot*cos(theta))/(3*cos(theta)^2 - 14);
end

