% System matrices
A = [1 2; 1 2];
B = [1; 1];

% Controllability matrix
Gamma = ctrb(A,B);

T = [[1; 1], [-1; 1]]

Ahat = inv(T)*A*T
Bhat = inv(T)*B