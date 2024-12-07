clear; close all
%% Kade Howes, ECE6320, Final Project, Problem 2 (Linear Control Design)
% A
data = load("sys2.mat");
A = data.A;
B = data.B;
C = data.C;
x_d = data.x_d;
n = height(A);

% check controllability
Gamma = ctrb(A,B);
rGamma = rank(Gamma); % not completely controllable
% check observability
Omega = obsv(A,C);
rOmega = rank(Omega); % not completely observable

% Kalman decomposition
% v1 = orth(Gamma);
% v2 = null(Gamma');
% v3 = null(Omega);
% v4 = orth(Omega');

% compute minimal realization
sys = ss(A,B,C,[]);
[msys, T] = minreal(sys);
Ahat = T*A*T^-1;
Bhat = T*B;
Chat = C*T^-1;
A11hat = Ahat(1:3,1:3);
n1 = height(A11hat);
B1hat = Bhat(1:3,:);
C1hat = Chat(:,1:3);

% check minimal realization
syms s;
TF1 = simplify(C*(s*eye(n)-A)^-1*B);
TF2 = simplify(C1hat*(s*eye(n1)-A11hat)^-1*B1hat);
diff = TF1-TF2;
diff_check1 = double(subs(diff,s,100)); % essentially 0
diff_check2 = double(subs(diff,s,103428)); % essentially 0
diff_checkA = A11hat-msys.A; % essentially 0
diff_checkB = B1hat-msys.B; % essentially 0
diff_checkC = C1hat-msys.C; % essentially 0

% convert and save names the problem asks for
Ahat = A11hat;
Bhat = B1hat;
Chat = C1hat;
save("prob2a.mat", "Ahat", "Bhat", "Chat");

% B
syms uff1 uff2 uff3 uff4;
uff = [uff1; uff2; uff3; uff4];
eqn = B*uff + A*x_d == 0; 
sol = solve(eqn, uff); % solve for uff
uff = double(subs(uff, uff, [sol.uff1; sol.uff2; sol.uff3; sol.uff4]));
% not controllable, do controllable decomp.
Tc = [orth(Gamma) null(Gamma')];
Achat = Tc^-1*A*Tc;
Bchat = Tc^-1*B; 
Ac22hat = Achat(rGamma+1:end,rGamma+1:end);
eigAc22hat = eig(Ac22hat); % stabilizable
Ac11hat = Achat(1:rGamma,1:rGamma); % controllable portion
Bc1hat = Bchat(1:rGamma,:);
K1 = place(Ac11hat, Bc1hat, -rGamma:-1);
K = [K1 zeros([width(B) n-rGamma])]*Tc^-1;
eigclosedc = eig(A-B*K);

% save results
save("prob2b.mat", "uff", "K");

% C
To = [orth(Omega') null(Omega)];
Aohat = To^-1*A*To;
Cohat = C*To;
Ao22hat = Aohat(rOmega+1:end,rOmega+1:end);
eigAo22hat = eig(Ao22hat); % not detectable, continue anyways
Ao11hat = Aohat(1:rOmega,1:rOmega); % observable portion
Co1hat = Cohat(:,1:rOmega);
L1 = (place(Ao11hat', Co1hat', -rOmega:-1))';
L = To*[L1; zeros([n-rOmega height(C)])];
eigclosedo = eig(A-L*C);

% save results
save("prob2c.mat", "L");