clear
close all;

% 2.3
g = 9.8; m = 1/9.8; l = 0.25; b = 1;
%a
A23a = [0 1; g/l (-b)/(m*l^2)];
eigA23a = eig(A23a);
B23a = [0; 1/(m*l^2)];
Gamma = [B23a A23a*B23a];
r_G23a = rank(Gamma); % =n => controllable and stabilizable

% b
A23b = [0 1; -g/l (-b)/(m*l^2)];
eigA23b = eig(A23b);
B23b = [0; 1/(m*l^2)];
Gamma = [B23b A23b*B23b];
r_G23b = rank(Gamma); % =n => controllable and stabilizable

% b
A23c = [0 1; -g/(l*sqrt(2)) (-b)/(m*l^2)];
eigA23c = eig(A23c);
B23c = [0; 1/(m*l^2)];
Gamma = [B23c A23c*B23c];
r_G23c = rank(Gamma); % =n => controllable and stabilizable

% 2.4
% c
k = 1;
A24 = [0 1; -k 0];
eigA24 = eig(A24);
B24 = [0; 1];
Gamma = [B24 A24*B24];
r_G24 = rank(Gamma); % =n => controllable and stabilizable

% 2.6
% b
g = 9.8; m = 1/9.8; I = 1; b = 2;
A26 = [0 1; (g*m)/I (-b)/I];
eigA26 = eig(A26);
B26 = [0; 1/I];
Gamma = [B26 A26*B26];
r_G26 = rank(Gamma); % =n => controllable and stabilizable

% 2.7
% b
A27 = [0 0 0; 0 0 0; 0 0 0];
eigA27 = eig(A27);
B27 = [1 0; 0 0; 0 1];
Gamma = [B27 A27*B27 A27^2*B27];
r_G27 = rank(Gamma); % =2 ~= n => not controllable
M26 = [eig(A27).*eye(3)-A27 B27];
r_M26 = rank(M26); % =2 ~= n => not stabilizable

% d
A27d = [0 1 0; -1 0 0; 0 0 0];
eigA27d = eig(A27d);
B27d = [1 -1; 0 0; 0 1];
Gamma = [B27d A27d*B27d A27d^2*B27d];
r_G27d = rank(Gamma); % =2 ~= n => not controllable
% M26d = [eig(A27d).*eye(3)-A27d B27d];
% r_M26d = rank(M26d); % =2 ~= n => not stabilizable

% Control Design
% simple systems
% S1
syms lam real
syms k1 k2
A1 = [1 2; 3 4]; B1 = [0; 1];
r_G1 = rank([B1 A1*B1]); % =2 => controllable
Abar1 = A1-(B1*[k1 k2]);
char1 = det(lam*eye(height(A1))-Abar1);
% check 
K = [5 7];
Abar1check = A1-B1*K;
eig1check = eig(Abar1check);
% S2
A2 = [-1 -2; 6 7]; B2 = [-0.5; 1];
r_G2 = rank([B2 A2*B2]); % =2 => controllable
Abar2 = A2-(B2*[k1 k2]);
char2 = det(lam*eye(height(A2))-Abar2);
eqn2 = [-0.5*k1+k2-6 == 2, 1.5*k1-2*k2+5 == 1];
sol2 = solve(eqn2);
K = [24 20];
Abar2check = A2-B2*K;
eig2check = eig(Abar2check);
% S3
syms k11 k12 k13 k21 k22 k23;
A3 = [1 2 3; 4 5 6; 7 8 9]; B3 = [1 0; 0 1; 1 1];
r_G3 = rank([B3 A3*B3 A3^2*B3]); % =3 => controllable
Abar3 = A3-B3*[k11 k12 k13; k21 k22 k23];
k11 = 3; k21 = 4;
k12 = 2; k22 = 6;
k13 = 0; k23 = 10;
K = [k11 k12 k13; k21 k22 k23];
Abar3check = A3-B3*K;
eig3check = eig(Abar3check);

% Orbit-plan motion
R = 200;
K = SatelliteControlDesign(R);


% Controllability proof

% Exact control design








