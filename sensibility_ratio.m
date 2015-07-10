
load_javaplex;
load 'C:\Users\Cássio\Dropbox\workspace\subspacestream\dados\D8\subspace_stream.mat';

vars = [10 15 20 25 30];

res = containers.Map('KeyType', 'int64', 'ValueType', 'any');

MSUP = 50;
MAXITERS = 2;
DIVISIONS = 1000;

tic
for v = vars
    res(v) = ptsstream(X, sup, 1000, 2, MSUP, MAXITERS, v, DIVISIONS);
end
toc

save sens_ratio.mat res 
