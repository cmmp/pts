
% p3c function
function assigns = p3c(X, sup)
    % step 1 - identify relevant sections in attributes
    S = findSections(X);
    % step 2 - find p-signatures with p = 2
    [P, A] = findPsignatures(X, S);
    % step 3 - make clusters:
    assigns = makeClusters(X, P, A);
end

function assigns = makeClusters(X, P, A)
    [N, ~] = size(X);
    assigns = zeros(N,1);
    MINPTS = 10;
    K = 0;
    
    for i = 1:size(A,1)
        ai = A(i,1); aj = A(i,2);
        pi = (X(:,ai) >= P(i,1)) & (X(:,ai) <= P(i,2));
        pj = (X(:,aj) >= P(i,3)) & (X(:,aj) <= P(i,4));
        pmix = find(pi & pj);
        % now use dbscan on pmix to separate any clusters that might
        % have been erroneously merged
        Xsub = X(pmix, [ai aj]);
        Eps = epsestimate(Xsub, MINPTS);
        [class, ~] = dbscan(Xsub, MINPTS, Eps);
        class(class == -1) = 0; % noise is 0 to me
        idx = find(class ~= -1);
        class(idx) = class(idx) + K;
        K = max(class);
        assigns(pmix) = class;
    end
end

% find p-signatures with p = 2
% X - data set 
% S - map from attributes to sections
% P - p-signatures matrix
% A - p-signatures attributes matrix
function [P, A] = findPsignatures(X, S)
    d = size(X, 2); %dimensionality
    P = zeros(d * d * 10, 4); % [lower_1 upper_1 lower_2 upper_2 ; ...]
    A = zeros(d * d * 10, 2); % [att_1 att_2 ; ...]
    nused = 0; % number of used combinations
    threshold = 1e-20; % threshold for poisson acceptance
    
    for i = 1:(d - 1)
        for j = (i+1):d
            maxJ = max(X(:,j)); minJ = min(X(:,j));
            totalWidth = maxJ - minJ;
            % test the combinations of sections from att_i and att_j:
            Si = S(i);
            Sj = S(j);
            Nsi = size(Si,1); Nsj = size(Sj,1);  
            for k = 1:Nsi
                Sik = Si(k,:);
                psik = (X(:,i) >= Sik(1)) & (X(:,i) <= Sik(2));
                SuppS = sum(psik); % support of S
                for l = 1:Nsj
                    % check the combination S_ik with Sjl
                    Sjl = Sj(l,:);
                    psjl =  (X(:,j) >= Sjl(1)) & (X(:,j) <= Sjl(2));
                    ObsS = sum(psik & psjl); % observed support of R = S U S'
                    width = (Sjl(2) - Sjl(1)) / totalWidth; % width of S'
                    ESupp = SuppS * width; % expected support of R = (S U S') / S
                    pd = makedist('Poisson', 'lambda', ESupp);
                    prob = pdf(pd, ObsS);
                    if ObsS > ESupp && prob < threshold
                        % found a cluster core
                        nused = nused + 1;
                        P(nused, 1:2) = Sik;
                        P(nused, 3:4) = Sjl;
                        A(nused,1) = i;
                        A(nused,2) = j;
                    end
                end
            end
        end
    end
    % remove extra spaces:
    P = P(1:nused,:);
    A = A(1:nused,:);
end

function S = findSections(X)
    d = size(X, 2);
    S = containers.Map('KeyType', 'int64', 'ValueType', 'any');
    
    for i = 1:d
        S(i) = analyzeAttribute(X(:,i));
    end
end

% analyze attribute y
function S = analyzeAttribute(y)
    [hc, edges] = histcounts(y, 'BinMethod', 'sturges');
    nb = length(hc); % number of bins
    S = zeros(nb, 2); % sections matrix --> [s_1l s_1u ; s_2l s2u ; ...]
    crit = chi2inv(1 - 0.001, nb); % critical value for chi2 dist.
    
    allUnif = false;
    
    a = ones(nb,1); % which bins are active
    nsec = 0; % number of used sections
    
    while ~allUnif
       % run chi-square test to check if hc is uniform
       
       hidx = find(a); % index of the subselections
       hcsub = hc(hidx); % sub selection
        
       msup = mean(hcsub);
       hstat = sum((hcsub - msup) .* (hcsub - msup) ./ msup);
       
       if hstat <= crit
           allUnif = true;
       else
           % pick the largest support bin, mark it and disable it
           nsec = nsec + 1;
           [M I] = max(hcsub);
           Ri = hidx(I); % real index
           a(Ri) = 0; % disable this bin
           % add it to sections:
           S(nsec,1) = edges(Ri); % lower value
           S(nsec,2) = edges(Ri + 1); % upper value          
       end
    end
    
    S = S(1:nsec,:); % remove unused space
    
    % we need to sort the sections considering their lower values:
    [Ss I] = sort(S(:,1));
    S = S(I,:);
    
    % now we need to merge adjacent bins:
    i = 1;
    while true
        if i >= size(S,1)
            break;
        end
        if S(i,2) == S(i+1,1) % upper this = lower next
            S(i,2) = S(i+1,2);
            S(i+1,:) = []; % delete this section
        else
            i = i + 1;
        end
    end
end
