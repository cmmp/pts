load_javaplex
load 'C:\Users\Cássio\workspace\subspacestream\dados\D8\subspace_stream.mat'
% load 'C:\Users\Cássio\Dropbox\workspace\subspacestream\dados\D16\subspace_stream.mat'
%load 'C:\Users\Cássio\Dropbox\workspace\subspacestream\dados\D32\subspace_stream.mat'
%  load 'C:\Users\Cássio\Dropbox\workspace\subspacestream\dados\D64\subspace_stream.mat'
% load 'C:\Users\Cássio\Dropbox\workspace\subspacestream\dados\D128\subspace_stream.mat'

MAXITERS = 2;
RATIO = 10;
DIVISIONS = 1000;

rng default;

Q = ptsstream(X, sup, 1000, 2, 50, MAXITERS, RATIO, DIVISIONS);
% save Q.mat Q

