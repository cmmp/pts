
load_javaplex;
load 'C:\Users\Cássio\Dropbox\workspace\subspacestream\dados\D8\subspace_stream.mat';

vars = [100 325 550 775 1000];

res = containers.Map('KeyType', 'int64', 'ValueType', 'any');

MSUP = 50;
MAXITERS = 2;
RATIO = 10;

tic
for v = vars
    res(v) = ptsstream(X, sup, 1000, 2, MSUP, MAXITERS, RATIO, v);
end
toc

save sens_divisions.mat res 
