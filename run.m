%load /home/cassio/work/workspace/subspacestream/subspace_stream.mat
load 'C:/Users/Cássio/work/workspace/subspacestream/subspace_stream.mat'
Xsub = X(2000:3000,:); supsub = sup(1:1000);
% plotmatrix(Xsub);
% return
tic
assigns = p3c(Xsub, supsub);
toc
gscatter(Xsub(:,1), Xsub(:,3), assigns);
