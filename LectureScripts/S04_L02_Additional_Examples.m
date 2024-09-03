% %% Problem 1
% % Define system matrices
% A = [3 6 4; 9 6 10; -7 -7 -9]
% B = 1/3*[-2 1; 1 -2; 1 1]
% C = [2 3 4; 2 1 3]
% 
% % Check for complete controllability
% Gamma = [B A*B A^2*B];
% rank_Gamma = rank(Gamma) % rank(Gamma) = 2 => not completely controllable
% 
% % Check for complete observability
% Omega = [C; C*A; C*A^2];
% rank_Omega = rank(Omega) % rank(Omega) = 3 => completely observable
% 
% % Create our controllability decomposition Transform
% basis_Gamma = orth(Gamma)
% compl_Gamma = null(Gamma')
% T = [basis_Gamma compl_Gamma]
% 
% % Construct our new matrices
% Abar = inv(T)*A*T
% Bbar = inv(T)*B
% Cbar = C*T
% 
% % Extract the minimal realization
% Amin = Abar(1:2, 1:2)
% Bmin = Bbar(1:2, :)
% Cmin = Cbar(:,1:2)
% 
% % Test realization
% syms s real
% G = C*inv(s*eye(3) - A)*B
% G_min = Cmin*inv(s*eye(2)-Amin)*Bmin
% 
% % Check difference
% diff = double(subs(G-G_min, s, 1000*rand - 500))

 
%% Problem 2
% Define system matrices
A = [3 6 4; 9 6 10; -7 -7 -9]
B = 1/3*[1 4; 4 1; -2 1]
C = [1 2 3; 3 3 6]

% Check for complete controllability
Gamma = [B A*B A^2*B];
rank_Gamma = rank(Gamma) % rank(Gamma) = 3 => completely controllable

% Check for complete observability
Omega = [C; C*A; C*A^2];
rank_Omega = rank(Omega) % rank(Omega) = 2 => not completely observable

% Create our controllability decomposition Transform
basis_Omega = null(Omega)
compl_Omega = orth(Omega')
T = [ compl_Omega basis_Omega]

% Construct our new matrices
Abar = inv(T)*A*T
Bbar = inv(T)*B
Cbar = C*T
 
% Extract the minimal realization
Amin = Abar(1:2, 1:2)
Bmin = Bbar(1:2, :)
Cmin = Cbar(:,1:2)

% Test realization
syms s real
G = C*inv(s*eye(3) - A)*B
G_min = Cmin*inv(s*eye(2)-Amin)*Bmin

% Check difference
diff = double(subs(G-G_min, s, 1000*rand - 500))



% %% Problem 3
% % Define system matrices
% A = [-1 0  0 0 0  0;
%       0 -2 1 0 0  0;
%       0 0 -1 0 0  0;
%       0 0 0 -3 0  0;
%       0 0 0  0 -3 1;
%       0 0 0  0 0 -1]
% B = [1 0; 0 0; 0 1; 0 1; 0 0; 1 0]
% C = [1 2 0 0 0 0;
%      0 0 0 1 1 0]
%  
%  % Check for complete controllability
% Gamma = ctrb(A,B);
% rank_Gamma = rank(Gamma) % rank(Gamma) = 5 => not completely controllable
% 
% % Check for complete observability
% Omega = obsv(A,C);
% rank_Omega = rank(Omega) % rank(Omega) = 4 => not completely observable
% 
% % Formulate a minimum realization
% sys = ss(A, B, C, []);
% [msys, T] = minreal(sys)
















