clear;
close all;

%% ECE 6320 HW4

%% 8.7
A7 = [0 1 0; 0 0 0; 0 0 -2];
n7 = height(A7);
syms s real
eAts7 = (s*eye(n7)-A7)^-1;
eAt7 = ilaplace(eAts7);

%% 8.8
A8 = [0 0 0; 0 -1 1; 0 0 -1];
n8 = height(A8);
syms s real
eAts8 = (s*eye(n8)-A8)^-1;
eAt8 = ilaplace(eAts8);

%% 8.13
A13 = [0 1; -1 -2];
n = height(A13);
B13 = [0; 1];
syms f1 real
syms f2 real
syms s real
f = [f1 f2];
char13 = det(s*eye(n)-(A13+B13*f));

f1 = 1;
f2 = 2;
f = [f1 f2];
D13 = eig(A13+B13*f);

%% 9.5
A5 = [-2 -1 0; 1 0 0; 0 1 0];
B5 = [1; 0; 0];
C5 = [0 1 0];
n = height(A5);
syms s real
G5 = C5*(s*eye(n)-A5)^-1*B5;
poles5 = roots(coeffs(G5^-1));

%% 9.6
A6 = [-2 0 0; 0 1 0; 0 0 -1];
B6 = [1; 0; -1];
C6 = [1 1 0];
D6 = 1;
n = height(A6);
syms s real

G6 = C6*(s*eye(n)-A6)^-1*B6 + D6;
poles6 = roots([1 2]); % (s+3)/(s+2)
D6 = eig(A6);

%% A
AA = [0 -2; 2 0];
BA = [0; 1];
n = height(AA);
% 1
eigA = eig(AA);
% 2
syms s real
GA = (s*eye(n)-AA)^-1*BA;
% 3
T = 0.1;
syms t real
AAbart = expm(AA*t);
AAbar = expm(AA*0.1);
BAbar = int(expm(AA*(T-t)),0,T) * BA;
BAbarcheck = -1*(0.5*[cos(2*0); sin(2*0)] - 0.5*[cos(2*0.1); sin(2*0.1)]);
% 4
eigAAbar = eig(AAbar);
stableCheckAAbar = abs(eigAAbar);
% 5
AAbarE = eye(n) + T*AA;
BAbarE = T*BA;
% 6
eigAAbarE = eig(AAbarE);
stableCheckAAbarE = abs(eigAAbarE);


%% B
AB = [-1 0; 0 -1];
BB = [0; 1];
n = height(AB);
% 1
eigB = eig(AB);
% 2
syms s real
GB = (s*eye(n)-AB)^-1*BB;
% 3
T = 0.1;
syms t real
ABbart = expm(AB*t);
ABbar = expm(AB*T);
BBbar = int(expm(AB*(T-t)),0,T) * BB;
BBbarcheck = [0; exp(0)]-[0; exp(-0.1)];
% 4
eigABbar = eig(ABbar);
stableCheckABbar = abs(eigABbar);
% 5
ABbarE = eye(n) + T*AB;
BBbarE = T*BB;
% 6
eigABbarE = eig(ABbarE);
stableCheckABbarE = abs(eigABbarE);

%% C
AC = [1 -2; 1 0];
BC = [1; 0];
n = height(AC);
% 1
eigC = eig(AC);
% 2
syms s real
GC = (s*eye(n)-AC)^-1*BC;
% 3
T = 0.1;
syms t real
ACbart = expm(AC*t);
ACbar = expm(AC*T);
BCbar = int(expm(AC*(T-t)),0,T) * BC;
% BCbarcheck = [0; exp(0)]-[0; exp(-0.1)];
% 4
eigACbar = eig(ACbar);
stableCheckACbar = abs(eigACbar);
% 5
ACbarE = eye(n) + T*AC;
BCbarE = T*BC;
% 6
eigACbarE = eig(ACbarE);
stableCheckACbarE = abs(eigACbarE);


%% D
AD = [-100 0; 0 -100];
BD = 0;
n = height(AD);
% 1
eigD = eig(AD);
% 2
syms s real
GD = (s*eye(n)-AD)^-1*BD;
% 3
T = 0.1;
syms t real
ADbart = expm(AD*t);
ADbar = expm(AD*T);
BDbar = int(expm(AD*(T-t)),0,T) * BD;
% BDbarcheck = [0; exp(0)]-[0; exp(-0.1)];
% 4
eigADbar = eig(ADbar);
stableCheckADbar = abs(eigADbar);
% 5
ADbarE = eye(n) + T*AD;
BDbarE = T*BD;
% 6
eigADbarE = eig(ADbarE);
stableCheckADbarE = abs(eigADbarE);

%% functions

