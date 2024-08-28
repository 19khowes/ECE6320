clear all
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