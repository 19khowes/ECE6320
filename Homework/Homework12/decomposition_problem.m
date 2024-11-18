function decomposition_problem()
    % Create system
    [A, B, C] = get_system()
    x_d = [-10/3; 3; -3; 43/12; 0]

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
    closeEigK = eig(A-B*K)

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


