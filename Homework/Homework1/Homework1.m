clear variables
close all
%% ECE6320 Homework 1, Kade Howes

% P1
P1 = [-3 5;7 -10] * [-1;3];
% P2
P2 = [4 5 1;3 7 10;1 0 1] * [1;2;3];
% P5
P5 = [1 4 2;0 0 0;1 0 9];
rankP5 = rank(P5);
% P6
A1 = [2 3 5;-4 2 3];
A2 = [1 0 1;5 2 1; 1 2 2];
A3 = [2 1 1; 1 1 0;1 0 1];
nullA1 = null(A1);
checkA1 = A1 * [-1/16; -13/8; 1];
nullA2 = null(A2);
nullA3 = null(A3);
checkA3 = A3 * [1; -1; -1];
% P7
colspA1 = colspace(sym(A1));
colspA2 = colspace(sym(A2));
colspA3 = colspace(sym(A3));

%P9
A2 = [1 0 1;5 2 1; 1 2 2];
[V_A2, D_A2] = eig(A2);
D_A2_hand = roots([-1 5 -5 10]);
A2_D2 = (A2 - eye(3)*D_A2_hand(2));
A2_D3 = (A2 - eye(3)*D_A2_hand(3));
A2_rank = rank(A2);

A3 = [2 1 1; 1 1 0;1 0 1];
[V_A3, D_A3] = eig(A3);
A3_rank = rank(A3);

A4 = [4 5 1; 3 7 10; 1 0 1];
[V_A4, D_A4] = eig(A4);
D_A4_hand = roots([-1 12 -23 56]);
A4_rank = rank(A4);

A5 = [4 0 0; 0 7 0; 0 0 1];
[V_A5, D_A5] = eig(A5);
A5_rank = rank(A5);

A6 = [2 1 1; 0 5 10; 0 0 8];
[V_A6, D_A6] = eig(A6);
A6_rank = rank(A6);

%% Unstable System Simulation
close all 
clear variables
Acalc = [2 0 0; 2 2 2; 3 0 -1];
[V_Ac, D_Ac] = eig(Acalc);

% initial state is eigenvector of magnitude 5 (negative eigenvalue, u=0)
SimulationScriptUnstable(5*V_Ac(:,2), @(t, x) 0);

% initial state is eigenvector of magnitude 5 (negative eigenvalue, u=1)
SimulationScriptUnstable(5*V_Ac(:,2), @(t, x) 1);

% initial state is eigenvector of magnitude 5 (positive eigenvalue, u=0)
SimulationScriptUnstable(5*V_Ac(:,1), @(t, x) 0);

% random initial state (u=0)
SimulationScriptUnstable(rand([3 1]), @(t, x) 0);

%% Nonlinear System Simulation
close all
clear variables
PendulumCart
