function Ex2_4
    P.A = [-0.038 18.984 0 -32.174; -0.001 -0.632 1 0; 0 -0.759 -0.518 0; 0 0 1 0];
    P.B = [10.1 0; 0 -0.0086; 0.025 -0.011; 0 0];
    x0 = [10; 0.1; 0.1; 0];
    % choose Q and R
    Q = diag(ones([1 height(P.A)])); 
    R = 1000*diag(ones([1 width(P.B)])); % prioritze small control
    P.K = lqr(P.A,P.B,Q,R);
    
    % simulate system
    dt = 0.1;
    t0 = 0;
    tf = 100;
    t_vec = t0:dt:tf;
    [tmat, xmat] = ode45(@(t,x)f(t,x,P), t_vec, x0);
    
    tmat = tmat'; xmat = xmat';
    umat = getControlVector(tmat, xmat, P);
    
    % plot system
    figure;
    sgtitle("Exercise 2.4");
    subplot(4,2,1); hold on;
    plot(t_vec, xmat(1,:));
    ylabel('V', 'fontsize', 12);
    title('States');
    
    subplot(4,2,3); hold on;
    plot(t_vec, xmat(2,:));
    ylabel('\alpha', 'fontsize', 12);
    xlabel('time');
    
    subplot(4,2,5); hold on;
    plot(t_vec, xmat(3,:));
    ylabel('q', 'fontsize', 12);
    
    subplot(4,2,7); hold on;
    plot(t_vec, xmat(4,:));
    ylabel('\theta', 'fontsize', 12);
    xlabel('time');
    
    subplot(4,2,4); hold on;
    plot(t_vec, umat(1,:));
    ylabel('\delta_{th}', 'fontsize', 12);
    title('Inputs');
    
    subplot(4,2,6); hold on;
    plot(t_vec, umat(2,:));
    ylabel('\delta_{e}', 'fontsize', 12);
    xlabel('time');
    
    
    function xdot = f(t,x,P) 
        u = -P.K*x; % should be x - xd??
        xdot = P.A*x + P.B*u;
    end
    
    function u_vec = getControlVector(tvec, xvec, P)
        len = size(tvec, 2);
        u_vec = zeros(2, len);
        for k = 1:len
            u_vec(:,k) = -P.K*xvec(:,k);
        end
    end
end

