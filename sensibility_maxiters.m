
load_javaplex;
load 'C:\Users\Cássio\Dropbox\workspace\subspacestream\dados\D8\subspace_stream.mat';

vars = 1:10;
%vars = [4];

res = containers.Map('KeyType', 'int64', 'ValueType', 'any');

MSUP = 50;
RATIO = 10;
DIVISIONS = 1000;

tic
for v = vars
    res(v) = ptsstream(X, sup, 1000, 2, MSUP, v, RATIO, DIVISIONS);
end
toc

save sens_maxiters.mat res 
