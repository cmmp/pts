load 'C:\Users\Cássio\workspace\subspacestream\dados\D8\subspace_stream.mat'

A = [[1     4];
     [2     7];
     [3     5];
     [3     7];
     [3     8];
     [5     8]
];

Xs = X(1:1000, [3 8]);
% plotmatrix(Xs)

histogram(Xs(:,2), 'binmethod', 'sturges')
hc = histcounts(Xs(:,2), 'binmethod', 'sturges');

PCRIT = 0.95;

xcrit = chi2inv(PCRIT, length(hc) - 1);

obs = sum((hc - mean(hc)) .^ 2);

fprintf(1, 'valor obs: %.4f, xcrit = %.4f\n', obs, xcrit);


 