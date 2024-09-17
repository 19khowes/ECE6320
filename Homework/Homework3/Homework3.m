clear;
close all;

%% ECE 6320 HW3

%% Realization and Equivalence

A1 = [1 2 3; 4 5 6; 7 8 9];
B1 = [1 0; 0 1; 1 1];
C1 = [1 1 0; 0 1 1];
G1 = transfer_function(A1, B1, C1);
lam1 = eig(A1);

A2 = [5 -6 -4; -8 9 7; -2 3 1];
B2 = [-1 -1; -1 0; 0 -1];
C2 = [0 1 -1; -1 0 -1];
G2 = transfer_function(A2, B2, C2);
lam2 = eig(A2);

A3 = [9 -7 8; -3 1 -2; 6 -4 5];
B3 = [-1 -1; 1 0; 0 -1];
C3 = [0 1 -1; -1 0 -1];
G3 = transfer_function(A3, B3, C3);
lam3 = eig(A3);

% verify algebraic equivalence
T = [0 1 0; 0 0 -1; -1 0 0];
Ti = T^-1;

A3check = Ti*A1*T;
B3check = Ti*B1;
C3check = C1*T;

%% Matrix Exponential and Cayley-Hamilton
syms t real
A1 = [1 1 0; 0 1 0; 0 0 1];
A1t = A1^t;
eA1t = expm(A1*t);

A2 = [1 1 0; 0 0 1; 0 0 1];
% A2t = A2^t;
eA2t = expm(A2*t);

A3 = [2 0 0 0; 2 2 0 0; 0 0 3 3; 0 0 0 3];
A3t = A3^t;
eA3t = expm(A3*t);

A31 = [2 0; 2 2];
A31t = A31^t;
A32 = [3 3; 0 3];
A32t = A32^t;

%% functions
function G = transfer_function(A, B, C)
    syms s real
    n = size(A,1);
    G = simplify(C*(s*eye(n) - A)^-1 * B);
end