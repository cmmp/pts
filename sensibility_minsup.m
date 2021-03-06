
load_javaplex;
load 'C:\Users\C�ssio\Dropbox\workspace\subspacestream\dados\D8\subspace_stream.mat';

msups = [10 20 30 40 50 60 70 80 90 100];
%msups = [50];

res = containers.Map('KeyType', 'int64', 'ValueType', 'any');

MAXITERS = 2;
RATIO = 10;
DIVISIONS = 1000;

for msup = msups
    res(msup) = ptsstream(X, sup, 1000, 2, msup, MAXITERS, RATIO, DIVISIONS);
end

save sens_msup.mat res 
