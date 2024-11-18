function LQRBryson
    A = [3 6 4; 9 6 10; -7 -7 -9]; B = [-2/3 1/3; 1/3 -2/3; 1/3 1/3];
    nA = height(A);
    Gamma = ctrb(A,B);
    v = orth(Gamma);
    w = null(Gamma');
    T = [v w];
    Ahat = T^-1*A*T;
    Bhat = T^-1*B; % not controllable or stabilizable

    % simulate
    % Create the time variables
    dt = 0.01
    t0 = 0;
    tf = 1;
    t_vec = t0:dt:tf;
    
    Q = diag([1/1^2 1/100^2 1/100^2]); 
    R = diag([1/5^2 1/10^2]);
    S = diag([1/10^2 1/1^2 1/2^2]);
    P_T = S;
    x0 = reshape(P_T,[],1);
    [tmat_r, xmat_r] = ode45(@(t,x)ricotti(t,x,A,B,Q,R), flip(t_vec), x0);
    xmat = flip(xmat_r);
    tmat = flip(tmat_r);
    P_mat_ode45 = reshape(xmat',nA,nA,[]);

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
        X = agg_mat(1:3, :);
        Y = agg_mat(4:6, :);
        
        % Calculate P(t)
        P_mat_exp(:,:, k) = Y/X;
    end

    % Plot the P values
    figure;
    sgtitle("LQR Bryson's Method: P elements")
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
    
    disp('')
    function x_dot = ricotti(t,x,A,B,Q,R)
        nA = height(A);
        R_inv = inv(R);
        P = reshape(x,nA,[]);
        P_dot = -A'*P - P*A - Q + P*B*R_inv*B'*P;
        x_dot = reshape(P_dot,[],1);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotPmat(t_vec, Pmat, color, linewidth, fontsize)
    [p11, p12, p13, p21, p22, p23, p31, p32, p33] = extractRows(Pmat);
    
    % Plot first row, first column
    subplot(3,3,1); hold on;
    plot(t_vec, p11, color, 'linewidth', linewidth);
    ylabel('p_{11}', 'fontsize', fontsize);
    
    % Plot first row, second column
    subplot(3,3,2); hold on;
    plot(t_vec, p12, color, 'linewidth', linewidth);
    ylabel('p_{12}', 'fontsize', fontsize);

    % Plot first row, second column (THIRD)
    subplot(3,3,2); hold on;
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

function [p11, p12, p13, p21, p22, p23, p31, p32, p33] = extractRows(P_mat)
    p11 = squeeze(P_mat(1,1,:));
    p12 = squeeze(P_mat(1,2,:));
    p13 = squeeze(P_mat(1,3,:));
    p21 = squeeze(P_mat(2,1,:));
    p22 = squeeze(P_mat(2,2,:));
    p23 = squeeze(P_mat(2,3,:));
    p31 = squeeze(P_mat(2,1,:));
    p32 = squeeze(P_mat(2,2,:));
    p33 = squeeze(P_mat(2,3,:));
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

