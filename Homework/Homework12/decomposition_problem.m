function decomposition_problem()
    % Create system
    [A, B, C] = get_system();
    x_d = [-10/3; 3; -3; 43/12; 0];

    % find uff
    syms u1 u2 u3;
    uff = [u1;u2;u3];
    eqn = B*uff + A*x_d == 0;
    uff_sol = solve(eqn,uff);
    uff = double([uff_sol.u1; uff_sol.u2; uff_sol.u3]);

    % controllability
    Gamma = ctrb(A,B);
    rGamma = rank(Gamma);
    v = orth(Gamma);
    w = null(Gamma');
    T = [v w];
    Ahat = T^-1*A*T;
    Bhat = T^-1*B;
    A11hat = Ahat(1:rGamma,1:rGamma);
    B1hat = Bhat(1:rGamma,:);
    % transform Q
    n = height(A); m = width(B);
    Q = diag(1:n); R = diag(1:m);
    Qhat = T'*Q*T;
    Q11hat = Qhat(1:rGamma,1:rGamma);
    K1hat = lqr(A11hat,B1hat,Q11hat,R);
    padding = zeros([m n-rGamma]);
    Khat = [K1hat padding];
    K = Khat*T^-1;
    closeEigK = eig(A-B*K);

    % observer design
    Omega = obsv(A,C);
    rOmega = rank(Omega);
    T = [orth(Omega') null(Omega)];
    Ahat = T^-1*A*T;
    Chat = C*T;
    A11hat = Ahat(1:rOmega,1:rOmega);
    C1hat = Chat(:,1:rOmega);
    L1hat = place(A11hat',C1hat',10*closeEigK(1:rOmega))';
    padding = zeros([n-rOmega, m]);
    Lhat = [L1hat; padding];
    L = T*Lhat;
    closeEigL = eig(A-L*C);

    P.A = A; P.B = B; P.C = C; P.K = K; P.L = L;
    P.xd = x_d; P.uff = uff;
    simulate(P);
end

function simulate(P)
    n = height(P.A);
    % initial conditions
    x0 = zeros([n 1]); % state
    xhat0 = ones([n 1]); % state estimate
    
    % Simulate forward in time
    t = 0:0.01:10;
    [tvec,xvec] = ode45(@(t,x)f(t,x,P), t, [x0; xhat0]);
    tvec = tvec'; xvec = xvec';
    
    % Plot 
    figure;
    sgtitle("Output Feedback w/ Decompositions");
    % Plot x1
    subplot(5,1,1);
    plot(tvec, xvec(n+1,:), 'k--', 'linewidth', 3); hold on;
    plot(tvec, xvec(1,:), 'linewidth', 2); 
    plot([tvec(1), tvec(end)], [P.xd(1) P.xd(1)], 'g:', 'linewidth', 2);
    xlabel('Time (s)');
    ylabel('x_1');
    legend('Estimate', 'State', 'Desired');

    % Plot x2
    subplot(5,1,2);
    plot(tvec, xvec(n+2,:), 'k--', 'linewidth', 3); hold on;
    plot(tvec, xvec(2,:), 'linewidth', 2); 
    plot([tvec(1), tvec(end)], [P.xd(2) P.xd(2)], 'g:', 'linewidth', 2);
    xlabel('Time (s)');
    ylabel('x_2');

    % Plot x3
    subplot(5,1,3);
    plot(tvec, xvec(n+3,:), 'k--', 'linewidth', 3); hold on;
    plot(tvec, xvec(3,:), 'linewidth', 2); 
    plot([tvec(1), tvec(end)], [P.xd(3) P.xd(3)], 'g:', 'linewidth', 2);
    xlabel('Time (s)');
    ylabel('x_3');

    % Plot x4
    subplot(5,1,4);
    plot(tvec, xvec(n+4,:), 'k--', 'linewidth', 3); hold on;
    plot(tvec, xvec(4,:), 'linewidth', 2); 
    plot([tvec(1), tvec(end)], [P.xd(4) P.xd(4)], 'g:', 'linewidth', 2);
    xlabel('Time (s)');
    ylabel('x_4');

    % Plot x5
    subplot(5,1,5);
    plot(tvec, xvec(n+5,:), 'k--', 'linewidth', 3); hold on;
    plot(tvec, xvec(5,:), 'linewidth', 2); 
    plot([tvec(1), tvec(end)], [P.xd(5) P.xd(5)], 'g:', 'linewidth', 2);
    xlabel('Time (s)');
    ylabel('x_5');

end

function xdot = f(t,x,P)
    % Extract states
    n = height(P.A);
    x_state = x(1:n);
    x_hat = x(n+1:end);

    % Calculate control
    u = -P.K*(x_hat-P.xd) + P.uff;

    % Calculate state dynamics
    x_dot_state = P.A*x_state + P.B*u;

    % Calculate measurement (y)
    y = P.C*x_state;

    % Calculate x_hat dynamics
    x_hat_dot = P.A*x_hat + P.B*u + P.L*(y-P.C*x_hat);

    % aggregate state dynamics
    xdot = [x_dot_state; x_hat_dot];
end

function [A, B, C] = get_system()
    A = ...
        [3.0000         0   -4.0000         0         0;
         0    2.2000    0.8000         0    3.6000;
         0         0   -1.0000         0         0;
    1.0000         0    4.0000    4.0000         0;
         0   -1.4000   -0.6000         0   -3.2000];

    B = ...
        [0    2.0000         0;
        0.8000    0.2000    0.6000;
        1.0000    1.0000         0;
             0   -1.0000         0;
       -0.6000   -0.4000   -0.2000];

    C = ...
     [1     0     0     1     0;
     0     0     1     1     0;
     0     2     0     0     1];
end


