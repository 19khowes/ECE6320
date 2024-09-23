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
zeros5 = roots(coeffs(G5^-1));


%% functions

