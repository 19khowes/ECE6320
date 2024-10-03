function prob2_3

% close all;

    A = [0, 1; -39.2, -156.8];

 

    % Choose a Q that is positive definite

    Q = eye(2); %[1, 2; 0, 3];

 

    %Solve this: A'*P + P*A = -Q

    P = lyap(A', Q)

    eig_P = eig(P)

    LE = A'*P + P*A

 

    % Calculate mu

    mu = -min(eig(Q))/max(eig(P))

 

    % Simulate forward in time

    delx_0 = [.1; .5];

    [tvec, xmat] = ode45(@(t,x) f_linear(t, x, A), [0 20], delx_0);

    xmat = xmat';

 

    % Plot

    v_t0 = delx_0'*P*delx_0;

    bound = 1/min(eig(P)) .* exp(mu.*tvec) .* v_t0;

    del_x_norm = zeros(1, length(tvec));

    for k = 1:length(tvec)

        del_x_norm(k) = norm(xmat(:,k))^2;

    end

 

    figure;

    plot(tvec, del_x_norm, 'b', 'LineWidth',3); hold on;

    plot(tvec, bound, 'r', 'LineWidth',3)

 

 

    % simulate the nonlinear system

    x_eq = [pi; 0];

    x_0 = delx_0 + x_eq;

    [tvec, xmat] = ode45(@(t,x) f_nonlinear(t,x), [0 20], x_0);

    xmat = xmat';

 

    % Convert back to delx

    xmat = xmat - x_eq;

 

 

    % Plot

    v_t0 = delx_0'*P*delx_0;

    bound = 1/min(eig(P)) .* exp(mu.*tvec) .* v_t0;

    del_x_norm = zeros(1, length(tvec));

    for k = 1:length(tvec)

        del_x_norm(k) = norm(xmat(:,k))^2;

    end

 

    figure;

    plot(tvec, del_x_norm, 'b', 'LineWidth',3); hold on;

    plot(tvec, bound, 'r', 'LineWidth',3)

 

   

end

 

 

function xdot = f_linear(t, x, A)

    xdot = A*x;

end

 

 

function xdot = f_nonlinear(t, x)

    g = 9.8;

    m = 1/9.8;

    l = 0.25;

    b = 1.;

 

    xdot = [x(2); g/l*sin(x(1))-b/(m*l^2)*x(2)];

end