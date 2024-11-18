function Modified2_1
%MODIFIED2_1 Summary of this function goes here
%   Detailed explanation goes here
    A = [0 1; -2 -3];
    B = [0; 1];
    Q = [1 0; 0 0]; R = 1; S = [1 0; 0 1];
    n = height(A);
    
    % Create the time variables
    dt = 0.1; % Increase this by 0.1 from 0.1 to 0.5 
    t0 = 0;
    tf = 10;
    t_vec = t0:dt:tf;
    
    % simulate using ode45
    P_T = S;
    x0 = reshape(P_T,[],1);
    [tmat_r, xmat_r] = ode45(@(t,x)ricotti(t,x,A,B,Q,R), flip(t_vec), x0);
    xmat = flip(xmat_r);
    tmat = flip(tmat_r);
    P_mat_ode45 = reshape(xmat',n,n,[]);
    
    % Exponential solution
    % Set initial conditions for X and Y
    X0 = eye(height(A));
    Y0 = S;
    
    % Create aggregate matrix
    M = [A -B*inv(R)*B'; -Q -A'];
    
    % Loop through and calculate P
    T = t_vec(end);
    for k = 1:length(t_vec)
        % Calculate X and Y
        t = t_vec(k);
        agg_mat = expm(M*(t-T)) * [X0; Y0];
        X = agg_mat(1:2, :);
        Y = agg_mat(3:4, :);
        
        % Calculate P(t)
        P_mat_exp(:,:, k) = Y/X;
    end
    
    % Plot the P values
    figure;
    sgtitle("Modified 2.1: P elements")
    % Plot the ode45 data
    linewidth = 2;
    color = 'b';
    fontsize = 12;
    plotPmat(t_vec, P_mat_ode45, color, linewidth, fontsize);
    
    % Plot the exponential data
    linewidth = 2;
    color = 'g:';
    plotPmat(t_vec, P_mat_exp, color, linewidth, fontsize);
    legend('ode45', 'numeric/exponential')
    
    % Plot the difference
    figure;
    sgtitle("Modified 2.1: P difference");
    P_diff = P_mat_ode45 - P_mat_exp;
    plotPmat(t_vec, P_diff, color, linewidth, fontsize);
    
    % Simulate system forward in time (ode45)
    x0 = [1; 2];
    [tmat45, xmat45] = ode45(@(t,x)f(t,x,A,B,Q,R,S,T), t_vec, x0);
    
    % Simulate system forward in time (euler)
    x = x0;
    xmatEuler = zeros(2,length(t_vec));
    xmatEuler(:,1) = x;
    for k = 1:length(t_vec)-1
        % get current time
        t = t_vec(k);
        % calculate derivative at time t
        x_dot = f(t,x,A,B,Q,R,S,T);
        % update next state
        x = x + x_dot*dt;
        % store state
        xmatEuler(:,k+1) = x;
    end
    
    % plot states
    figure;
    sgtitle("Modified 2.1: States over time");
    color = 'b';
    plotStates(tmat45', xmat45', color, linewidth, fontsize);
    color = 'r:';
    plotStates(tmat45', xmatEuler, color, linewidth, fontsize);
    legend('ode45', 'Euler');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Dynamic functions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function x_dot = f(t,x,A,B,Q,R,S,T)
        P = calculateP(t,T,A,B,Q,R,S);
        u = -inv(R)*B'*P*x;
        x_dot = A*x + B*u;
    end
    
    function x_dot = ricotti(t,x,A,B,Q,R)
        n = height(A);
        R_inv = inv(R);
        P = reshape(x,n,[]);
        P_dot = -A'*P - P*A - Q + P*B*R_inv*B'*P;
        x_dot = reshape(P_dot,[],1);
    end
    
    function P = calculateP(t,T,A,B,Q,R,S)     
        % Set initial conditions for X and Y
        X0 = eye(2);
        Y0 = S;
        
        % Create aggregate matrix
        M = [A -B*R*B'; -Q -A'];
        
        agg_mat = expm(M*(t-T)) * [X0; Y0];
        X = agg_mat(1:2, :);
        Y = agg_mat(3:4, :);
        
        % Calculate P
        P = Y/X;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plotting functions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function plotPmat(t_vec, Pmat, color, linewidth, fontsize)
        [p11, p12, p21, p22] = extractRows(Pmat);
        
        % Plot first row, first column
        subplot(4,1,1); hold on;
        plot(t_vec, p11, color, 'linewidth', linewidth);
        ylabel('p_{11}', 'fontsize', fontsize);
        
        
        % Plot first row, second column
        subplot(4,1,2); hold on;
        plot(t_vec, p12, color, 'linewidth', linewidth);
        ylabel('p_{12}', 'fontsize', fontsize);
        
        % Plot second row, first column
        subplot(4,1,3); hold on;
        plot(t_vec, p21, color, 'linewidth', linewidth);
        ylabel('p_{21}', 'fontsize', fontsize);
        
        % Plot second row, second column
        subplot(4,1,4); hold on;
        plot(t_vec, p22, color, 'linewidth', linewidth);
        ylabel('p_{22}', 'fontsize', fontsize);
        xlabel('time');
    end
    
    function [p11, p12, p21, p22] = extractRows(P_mat)
        p11 = squeeze(P_mat(1,1,:));
        p12 = squeeze(P_mat(1,2,:));
        p21 = squeeze(P_mat(2,1,:));
        p22 = squeeze(P_mat(2,2,:));
    end
    
    function plotStates(t_vec, x_mat, color, linewidth, fontsize)
        % Plot first row, first column
        subplot(2,1,1); hold on;
        plot(t_vec, x_mat(1,:), color, 'linewidth', linewidth);
        ylabel('x_{1}', 'fontsize', fontsize);
        
        % Plot first row, second column
        subplot(2,1,2); hold on;
        plot(t_vec, x_mat(2,:), color, 'linewidth', linewidth);
        ylabel('x_{2}', 'fontsize', fontsize);
        xlabel('time');
    end
end

