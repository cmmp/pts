load_javaplex

javaaddpath('C:\Users\C�ssio\workspace\gridseg\target\gridseg-0.0.1-SNAPSHOT-jar-with-dependencies.jar');

load 'C:\Users\C�ssio\workspace\subspacestream\dados\D8\subspace_stream.mat'
% load 'C:\Users\C�ssio\Dropbox\workspace\subspacestream\dados\D16\subspace_stream.mat'
%load 'C:\Users\C�ssio\Dropbox\workspace\subspacestream\dados\D32\subspace_stream.mat'
%  load 'C:\Users\C�ssio\Dropbox\workspace\subspacestream\dados\D64\subspace_stream.mat'
% load 'C:\Users\C�ssio\Dropbox\workspace\subspacestream\dados\D128\subspace_stream.mat'

MAXITERS = 2;
RATIO = 10;
DIVISIONS = 1000;

rng default;

Q = ptsstream(X, sup, 1000, 2, 50, MAXITERS, RATIO, DIVISIONS);
% save Q.mat Q

