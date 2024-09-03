function S05_L03_in_class
    %K = controlDesign()
    %K = integralControlDesign()
    K = doubleIntegralControlDesign();

    x0 = [3;2;1];
    u = @(x) -K*x;

    % Simulate forward
%     [tmat, xmat] = ode45(@(t,x)f(t, x, u(x)), [0 10], x0);
%     xmat = xmat';

%     [tmat, xmat] = ode45(@(t,z)f_integral(t, z, u(z)), [0 30], [x0; 0;0]);
%     xmat = xmat';

    [tmat, xmat] = ode45(@(t,z)f_double_integral(t, z, u(z)), [0 30], [x0; 0;0; 0]);
    xmat = xmat';


    % Plot the states over time
    figure;
    subplot(4,1,1);
    plot(tmat, xmat(1,:));

    subplot(4,1,2);
    plot(tmat, xmat(2,:));

    subplot(4,1,3);
    plot(tmat, xmat(3,:));

    subplot(4,1,4);
    plot(tmat, xmat(4,:));
end

function K = controlDesign()
    A = [0 0 1; 1 2 0; 0 0 0];
    B = [0 1; 0 0 ; 1 0];

    Q = diag([1, 3, 0]);
    R = diag([2, 3]);

    P = are(A, B*inv(R)*B', Q);
    K = lqr(A, B, Q, R);
end

function K = integralControlDesign()
    A = [0 0 1; 1 2 0; 0 0 0];
    B = [0 1; 0 0 ; 1 0];

    % Form our augmented system
    G = [1 0 0; 0 0 1];
    Abar = [A zeros(3,2); G zeros(2,2)];
    Bbar = [B; zeros(2, 2)];

    % Evaluate the controllability
    assert(rank(ctrb(Abar, Bbar)) == 5, "not controllable");

    % Design our Q and R matrices
    Q = diag([1, 3, 0, .1, .1]);
    R = diag([2, 3]);

    % Evaluate the observability condition
    Cc = sqrt(Q);
    
    assert(rank(obsv(Abar, Cc))== 5, "Not observable");


    K = lqr(Abar, Bbar, Q, R);
end

function K = doubleIntegralControlDesign()
    A = [0 0 1; 1 2 0; 0 0 0];
    B = [0 1; 0 0 ; 1 0];

    % Form our augmented system
    G = [1 0 0; 0 0 1];
    G2 = [0 1];
    %          sigma       gamma     x   sigma     gamma          x       sigma  gamma
    Abar = [A zeros(3,2) zeros(3,1); G zeros(2,2) zeros(2,1); zeros(1,3), G2, zeros(1,1)];
    Bbar = [B; zeros(2, 2); zeros(1,2)];

    % Evaluate the controllability
    assert(rank(ctrb(Abar, Bbar)) == 6, "not controllable");

    % Design our Q and R matrices
    Q = diag([1, 3, 0, 10, .1, .1]);
    R = diag([2, 3]);

    % Evaluate the observability condition
    Cc = sqrt(Q);
    
    assert(rank(obsv(Abar, Cc))== 6, "Not observable");


    K = lqr(Abar, Bbar, Q, R);
end

function xdot = f(t, x, u)
    A = [0 0 1; 1 2 0; 0 0 0];
    B = [0 1; 0 0 ; 1 0];

    w = [9; 0; sin(5*t)];

    xdot = A*x + B*u + w;
end

function zdot = f_integral(t, z, u)
    xdot = f(t, z(1:3), u);

    G = [1 0 0; 0 0 1];
    sigma_dot = G*z(1:3); % reality, you keep track of this based on x


    zdot = [xdot; sigma_dot];    
end

function zdot = f_double_integral(t, z, u)
    xdot = f(t, z(1:3), u);

    G = [1 0 0; 0 0 1];
    sigma_dot = G*z(1:3); % reality, you keep track of this based on x
    gamma_dot = z(5); % sigma 2

    zdot = [xdot; sigma_dot; gamma_dot];    
end