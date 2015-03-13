%load /home/cassio/work/workspace/subspacestream/subspace_stream.mat
load 'C:/Users/Cássio/work/workspace/subspacestream/subspace_stream.mat'
Xsub = X(2000:3000,:); supsub = sup(1:1000);
% plotmatrix(Xsub);
% return
tic
assigns = p3c(Xsub, supsub);
toc
gscatter(Xsub(:,1), Xsub(:,4), assigns);

% cluster 1 -> problematic in (1,4)

% Xsub = Xsub(find(assigns == 1),[1 4]);
% % plotmatrix(Xsub)
% x = Xsub;
% [m,n]=size(x);
% k = 10;
% Eps=((prod(max(x)-min(x))*k*gamma(.5*n+1))/(m*sqrt(pi.^n))).^(1/n);
% fprintf('eps found: %.4f\n', Eps);
% [class, type] = dbscan(x, k, Eps);
% gscatter(x(:,1), x(:,2), class);


% Xsub = Xsub(find(assigns == 2),[2 4]);
% % plotmatrix(Xsub)
% x = Xsub;
% [m,n]=size(x);
% k = 10;
% Eps=((prod(max(x)-min(x))*k*gamma(.5*n+1))/(m*sqrt(pi.^n))).^(1/n);
% fprintf('eps found: %.4f\n', Eps);
% [class, type] = dbscan(x, k, Eps);
% gscatter(x(:,1), x(:,2), class);
