function Q = ptsstream(X, sup, wsize, K, minsuprt, maxiters, ratio, divisions)
% Streaming algorithm
% Q = ptsstream(X, sup)
%
% PARAMETERS:
% X - data matrix (m,n) with m instances of n observed variables
% sup - array of m elements containing numeric supervision
% wsize - window size
% K - how many final topological clusters are desired
% minsuprt - minimum number of points to consider a psignature relevant.
% maxiters - maximum number of EM iterations
% ratio - ratio for landmark selection
% divisions - number of divisions for defining the step of epsilon
%
% RETURN:
% Q - array of m elements as measured according to ARI evaluation
%     measure.
%
% NOTES:
%  Current implementation only works for finding 2-d embedded clusters.
%
    [N, d] = size(X);
    
    if mod(N, wsize) > 0
        error('please use a wsize that divides the data set size...');
    end        
        
    Q = zeros(N / wsize, 1); % quality measurements
    
    j = 1;
    for i = 1:wsize:N
        slice = i:(i+wsize-1);
        tic
        Q(j) = ptsrun(X(slice,:), sup(slice), 'on', K, minsuprt, maxiters, ratio, divisions);
        toc
        memory
        fprintf('Got ARI = %.2f for window [%d:%d/%d]\n', ...
            Q(j), i, i + wsize - 1, N);
        j = j + 1;
    end 
end