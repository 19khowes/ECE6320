%% Example 1:
V = [1 0 1; 0 1 0; -2 0 -2];
%V = [1 4 2; 0 0 0; 1 0 9]
N_V = null(V) % Compute the null space
im_V = orth(V) % Compute the image space
[vectors,values] = eig(V) % Compute the eigenvectors and eigenvalues
eig_values = diag(values)
