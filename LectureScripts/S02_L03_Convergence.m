function S02_L03_Convergence
    % Initialize variables
    x0 = [-2; -2];
    t = 0:.01:20;
    
    % Define the matrices
    A = [0 1; -1 -2];
    Q = diag([1, 1]);
    P = lyap(A', Q);
    
    % Get the convergence rate
    lam_Q_min = min(eig(Q));
    lam_P_max = max(eig(P));
    lam_P_min = min(eig(P));
    rate = lam_Q_min / lam_P_max
    v0 = x0'*P*x0;
    
    % Simulate forward in time
    [tvec xvec] = ode45(@(t, x) f(t, x, A), t, x0);
    xvec = xvec';
    
    % Calculate norm and upper bound
    norm_vec = zeros(size(t));
    upper_vec = zeros(size(t));
    for k = 1:length(t)
        norm_vec(k) = norm(xvec(:,k))^2;
        upper_vec(k) = 1/lam_P_min * exp(-rate*tvec(k))*v0;
    end
    
    % Plot the result
    figure;
    plot(tvec, norm_vec, 'b', 'linewidth', 3); hold on;
    plot(tvec, upper_vec, 'r', 'linewidth', 3);
    legend('norm', 'Upper bound');
    set(gca, 'fontsize', 12);    
end

function xdot = f(t,x, A)    
    xdot = A*x;
end