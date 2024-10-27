function control_decomposition()
    control_design_system_1()
    control_design_system_2()
    control_design_system_3()

end

function control_design_system_1()
    % Get system
    [A, B] = get_system1();
    display("Create your controller for system 1");
    n = height(A);
    m = width(B);
    Gamma = ctrb(A,B);
    rGamma = rank(Gamma); % rank of Gamma is 3 ~= 5 (not fully controllable)
    v = orth(Gamma);
    w = null(Gamma');
    T = [v w];
    Ahat = T^-1*A*T;
    Bhat = T^-1*B;
    controlDim = 3; % size of Ahat_11 (number of placeable eigenvalues)
    Q = diag(1:controlDim);
    R = diag(1:m);
    Ac = Ahat(1:controlDim,1:controlDim); % controllable portion of Ahat
    Bc = Bhat(1:controlDim,:); % controllable portion of Bhat
    Khat = lqr(Ac, Bc, Q, R); % control for controllable portion
    Khat = [Khat zeros([m n-controlDim])];
    K = Khat*T^-1
    Abar = A-B*K;
    eigAbar = eig(Abar)
end

function control_design_system_2()
    % Get system
    [A, B] = get_system2();
    display("Create your controller for system 2");    
end

function control_design_system_3()
    % Get system
    [A, B] = get_system3();
    display("Create your controller for system 3");    
end


function [A, B] = get_system1()
    A = ...
       [-2.200         0   -0.4000         0         0;
             0    1.0000         0   -1.0000    1.0000;
        0.6000         0   -0.8000         0         0;
             0    1.0000         0    3.0000         0;
             0    3.0000         0    2.0000    2.0000];

    B = ...
        [0     0;
         1     0;
         0     0;
        -1     1;
         0     0];
end

function [A, B] = get_system2()

    A = ...
        [2.2000         0    0.4000         0         0;
         0    2.0000         0   -1.0000    3.0000;
   -0.6000         0    0.8000         0         0;
         0    2.0000         0    5.0000   -2.0000;
         0    2.0000         0    2.0000    2.0000];

    B = ...
        [0         0   -0.2000;
         0   -1.0000         0;
         0         0    0.6000;
         0    2.0000         0;
    1.0000    1.0000         0];
end

function [A, B] = get_system3()
    A = ...
        [5.0000         0   -2.0000    2.0000         0;
        0.4000    1.6000    0.8000    0.4000    1.8000;
        2.0000         0    2.0000    2.0000         0;
       -1.0000         0    3.0000    2.0000         0;
       -0.8000   -1.2000   -1.6000   -0.8000   -2.6000];

    B = ...
        [0    2.0000         0;
        0.8000    0.2000    0.6000;
        1.0000    1.0000         0;
             0   -1.0000         0;
       -0.6000   -0.4000   -0.2000];
end

