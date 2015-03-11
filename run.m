%load /home/cassio/work/workspace/subspacestream/subspace_stream.mat
load 'C:/Users/Cássio/work/workspace/subspacestream/subspace_stream.mat'
Xsub = X(1:1000,:); supsub = sup(1:1000);
tic
p3c(Xsub, supsub);
toc