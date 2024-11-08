%% ECE6320, HW10, Smooth Differential Drive
clear; close all;
x0 = [0; 0; 0; 0; 0];

% Define parameters
P.r = 1;
P.L = 2;
P.ep = 1;
P.A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0];
P.B = [0 0; 0 0; 1 0; 0 1];
P.K = place(P.A,P.B,[-1 -2 -3 -4]);

% Simulate the state forward in time the state

dt = 0.01;
t = [0:dt:10];
% xd = xeq;
% u = @(t,x,P) ueq-K*(x-xeq);
% u = @(t,x,P) [1;2];
u = @(t,x,P) u_function(t,x,P);

[tmat, xmat] = ode45(@(t,x)f(t,x,u,P), t, x0);
tmat = tmat';
xmat = xmat';
umat = getControlVector(tmat, xmat, u, P);

y_ep = xmat(1:2,:) + P.ep*[cos(xmat(3,:)); sin(xmat(3,:))];
qd = [sin(t); t];

%% Plot the results
fontsize = 12;

% Plot the states
figure;
subplot(5,1,1);
plot(tmat, xmat(1,:), 'b', 'linewidth', 3); hold on;
plot(tmat, y_ep(1,:), '--r', 'linewidth', 2);
plot(tmat, qd(1,:), '--g', 'linewidth', 2);
ylabel('$x_1(t)$, $x_{\epsilon 1}(t)$, $q_{d1}(t)$', 'fontsize', fontsize, 'Interpreter','latex');
set(gca, 'fontsize', fontsize);
legend('x_1', 'x_{\epsilon 1}', 'q_{d1}')

subplot(5,1,2);
plot(tmat, xmat(2,:), 'b', 'linewidth', 3); hold on;
plot(tmat, y_ep(2,:), '--r', 'linewidth', 2);
plot(tmat, qd(2,:), '--g', 'linewidth', 2);
ylabel('$x_2(t)$, $x_{\epsilon 2}(t)$, $q_{d2}(t)$', 'fontsize', fontsize, 'Interpreter','latex');
set(gca, 'fontsize', fontsize);
legend('x_2', 'x_{\epsilon 2}', 'q_{d2}')

subplot(5,1,3);
plot(tmat, xmat(3,:), 'b', 'linewidth', 3);
ylabel('$\theta(t)$', 'fontsize', fontsize, 'Interpreter','latex');
set(gca, 'fontsize', fontsize);

subplot(5,1,4);
plot(tmat, xmat(4,:), 'b', 'linewidth', 3);
ylabel('$w_r(t)$', 'fontsize', fontsize, 'Interpreter','latex');
set(gca, 'fontsize', fontsize);

subplot(5,1,5);
plot(tmat, xmat(5,:), 'b', 'linewidth', 3);
ylabel('$w_l(t)$', 'fontsize', fontsize, 'Interpreter','latex');
set(gca, 'fontsize', fontsize);

% Plot the input
figure;
subplot(2,1,1);
plot(tmat, umat(1,:), 'b', 'linewidth', 3);
ylabel('$\dot{w_r}(t)$', 'fontsize', fontsize, 'Interpreter','latex');
set(gca, 'fontsize', fontsize);

subplot(2,1,2);
plot(tmat, umat(2,:), 'b', 'linewidth', 3);
ylabel('$\dot{w_l}(t)$', 'fontsize', fontsize, 'Interpreter','latex');
set(gca, 'fontsize', fontsize);

%% input
function abar = u_function(t,x,P)
    r = P.r;
    L = P.L;
    % get states
    q1 = x(1);
    q2 = x(2);
    theta = x(3);
    wr = x(4);
    wl = x(5);

    v = (r/2)*(wr+wl);
    w = (r/L)*(wr-wl);

    M = [r/2 r/2; r/L -r/L];
    Re = [cos(theta) -P.ep*sin(theta); sin(theta) P.ep*cos(theta)];
    what = [0 -P.ep*w; w/P.ep 0];
    vbar = [v; w];

    z = [q1+P.ep*cos(theta);
         q2+P.ep*sin(theta);
         v*cos(theta)-P.ep*w*sin(theta);
         v*sin(theta)+P.ep*w*cos(theta);];

    zd = [sin(t); t; cos(t); 1];
    
    u = [-sin(t); 0] - P.K*(z-zd);

    abar = Re^-1*u - what*vbar;
    abar = M^-1 * abar;
end

function u_vec = getControlVector(tvec, xvec, u, P)
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
        u_vec(:,k) = u(tvec(k), xvec(:,k), P);
    end

end

%% dynamics
function xdot = f(t,x,u_function,P)
    r = P.r;
    L = P.L;
    
    % get states
    theta = x(3);
    wr = x(4);
    wl = x(5);
    
    v = (r/2)*(wr+wl);
    w = (r/L)*(wr-wl);

    % get inputs
    u = u_function(t,x,P);

    xdot = [v*cos(theta); v*sin(theta); w; u(1); u(2)];
end

