load_javaplex
% load 'C:\Users\Cássio\Dropbox\workspace\subspacestream\dados\D8\subspace_stream.mat'
load 'C:\Users\Cássio\Dropbox\workspace\subspacestream\dados\D16\subspace_stream.mat'
% load 'C:\Users\Cássio\Dropbox\workspace\subspacestream\dados\D32\subspace_stream.mat'
% load 'C:\Users\Cássio\Dropbox\workspace\subspacestream\dados\D64\subspace_stream.mat'
% load 'C:\Users\Cássio\Dropbox\workspace\subspacestream\dados\D128\subspace_stream.mat'

Q = p3cstream(X, sup, 1000, 2, 50);
% save Q.mat Q

