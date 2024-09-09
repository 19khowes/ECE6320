clear variables
close all
%% ECE6320 Homework 2, Kade Howes

% 2.3
% a
g = 9.8; m = 1/9.8; l = 0.25; b = 1;
Aa = [0 1; g/l (-b)/(m*l^2)];
ea = eig(Aa);
eacheck = roots([1 156.8 -39.2]);
% b
g = 9.8; m = 1/9.8; l = 0.25; b = 1;
Ab = [0 1; -g/l (-b)/(m*l^2)];
eb = eig(Ab);
ebcheck = roots([1 156.8 39.2]);
% b
g = 9.8; m = 1/9.8; l = 0.25; b = 1;
Ac = [0 1; g/(l*sqrt(2)) (-b)/(m*l^2)];
ec = eig(Ac);
eccheck = roots([1 156.8 - ...
    Ac(2,1)]);

