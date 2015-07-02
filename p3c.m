
function J = p3c(X, sup, topo, K, minsuprt)
load_javaplex;
% p3c-based algorithm
% CALL: assigns = p3c(X, sup)
%
% PARAMETERS:
% X - (m,n) with m instances and n variables input matrix
% sup - (m) array with supervision
% topo - on/off whether to run topo analysis
% K - how many topological clusters are desired
%
% RETURNS: 
% J - jaccard index value

    % step 1 - identify relevant sections in attributes
    S = findSections(X);
    % step 2 - find p-signatures with p = 2
    [P, A] = findPsignatures(X, S, minsuprt);
    % step 3 - make clusters:
%     [assigns, A] = makeClusters(X, P, A);
    assigns = makeClustersEM(X, P, A, minsuprt);
    
%     % step 4 - run topological analysis?
    if strcmp(topo, 'on')
        assigns = topoAnalysis(X, A, assigns, K);
    end
    % compute jaccard:
    %J = jaccard(assigns, sup);
    J = MyUtils.computeAdjustedRandIndex(int32(assigns), int32(sup));
%     clf(); plotWindow(X, assigns);
end

function plotWindow(X, assigns)
%    set(0,'DefaultTextFontSize',26);
    ngroups = length(unique(assigns));
    %cmap = lines(ngroups);
    %cmap = myrandomjetcmap();
    cmap = [];
    PSIZE = 12;
    subplot(2,3,1);
    Xsub = X;
    gscatter(Xsub(:,1),Xsub(:,2), assigns, cmap, [], PSIZE, 'on'); xlabel('1'); ylabel('2');
    subplot(2,3,2);
    gscatter(Xsub(:,1),Xsub(:,3), assigns, cmap, [], PSIZE, 'on'); xlabel('1'); ylabel('3');
    subplot(2,3,3);
    gscatter(Xsub(:,1),Xsub(:,4), assigns, cmap, [], PSIZE, 'on'); xlabel('1'); ylabel('4');
    subplot(2,3,4);
    gscatter(Xsub(:,2),Xsub(:,3), assigns, cmap, [], PSIZE, 'on'); xlabel('2'); ylabel('3');
    subplot(2,3,5);
    gscatter(Xsub(:,2),Xsub(:,4), assigns, cmap, [], PSIZE, 'on'); xlabel('2'); ylabel('4');
    subplot(2,3,6);
    gscatter(Xsub(:,3),Xsub(:,4), assigns, cmap, [], PSIZE, 'on'); xlabel('3'); ylabel('4');
end

function topoassigns = topoAnalysis(X, A, assigns, K)
    
	% Constant for topological analysis, maximum dimension */
	MAX_D = 3;
	% Ratio for selection of landmarks -> 1 landmark : ratio points */
	RATIO = 10;
	% Number of divisions for filtration analysis */
	DIVISIONS = 1000;

    k = size(A,1); % number of actual clusters
    KMX = zeros(k, MAX_D); % kmeans input matrix
    
    topoassigns = zeros(size(assigns));
    
    for i = 1:k
        pts = find(assigns == i);
        npts = length(pts);
        if npts <= 1
            KMX(i,:) = -10;
            continue;
        end
        
        Xproj = X(pts, [A(i,1) A(i,2)]);
        
        maxdist = max(pdist(Xproj));
        ems = edu.stanford.math.plex4.metric.impl.EuclideanMetricSpace(Xproj);
        
        if npts < RATIO
            nLandmarks = npts;
        else
            nLandmarks = round(npts / RATIO);
        end
        
        maxmin = edu.stanford.math.plex4.metric.landmark.MaxMinLandmarkSelector(ems, nLandmarks);
        lt = edu.stanford.math.plex4.streams.impl.LazyWitnessStream(ems, maxmin, MAX_D, maxdist, 0, DIVISIONS);
        lt.finalizeStream();
        persistence = edu.stanford.math.plex4.api.Plex4.getModularSimplicialAlgorithm(MAX_D, 2);
        bc = persistence.computeIntervals(lt);
        nd = length(bc.getBettiSequence());
        KMX(i,1:nd) = bc.getBettiSequence();
    end
    
    if any(KMX(:,1) == -10) % we have non-significant clusters
        K = K + 1; 
    end
    
    KMX = KMX(:,2:end); % remove 0-d points
    tclusters = kmeans(KMX, K); % topological clusters
    
    for i = 1:k
        topoassigns(assigns == i) = tclusters(i);
    end
end

function [assigns] = makeClustersEM(X, P, A, minsuprt)
% this function creates the clusters based on the EM algorithm
% described in the P3C paper along with some mods by us.
    
    n = size(X, 1);
    k = size(A, 1);
    M = zeros(n, k); % fuzzy membership matrix
    
    realassigns = zeros(n,1);
    
    S = containers.Map('KeyType', 'int64', 'ValueType', 'any'); %  covariances matrices map
    SI = containers.Map('KeyType', 'int64', 'ValueType', 'any'); %  inverse covariances matrices map
    C = zeros(k,2); % centroids matrix
 
    % compute basic covariances and means for each cluster:
    for i = 1:k
        ai = A(i,1); aj = A(i,2);
        pi = (X(:,ai) >= P(i,1)) & (X(:,ai) <= P(i,2));
        pj = (X(:,aj) >= P(i,3)) & (X(:,aj) <= P(i,4));
        pmix = pi & pj;
        M(:,i) = pmix;
        Xsub = X(pmix, [ai aj]);
        
        S(i) = cov(Xsub);
        SI(i) = inv(S(i));
        C(i,:) = mean(Xsub);
    end
    
    % remove the points that do not belong to any p-signature.
    % we just regard them as noise.
    
    Msum = sum(M,2);
    
    noise = Msum == 0;
    
    X = X(~noise,:);
    M = M(~noise,:);
    n = size(X, 1);
    
    % initialize matrix M by assigning each point to the 'closest'
    % cluster core in terms of mahalanobis distances to means of support
    % sets of cluster cores
    M(:) = 0;
    
    for i = 1:n
        minDist = inf; bestJ = -1;
        for j = 1:k
            ai = A(j,1); aj = A(j,2);
            x = X(i, [ai aj])';
            mu = C(j,:)';
            d = sqrt((x-mu)' * SI(j) * (x-mu)); % mahalanobis distance
            if d < minDist
                minDist = d;
                bestJ = j;
            end
        end
        M(i,bestJ) = 1;
    end
    
    % set initial alpha values:
    alphas = sum(M, 1) ./ n;  % mixture weights --> priors of a point coming from the i-th gaussian
    
    % recompute initial centroids and covariances:
    for i = 1:k
        ai = A(i,1); aj = A(i,2);
        Xsub = X(logical(M(:,i)),[ai aj]);
        if(size(Xsub,1) == 0)
            %error('got no points for this cluster!!');
             continue;
        end
        C(i,:) = mean(Xsub);
        S(i) = cov(Xsub); 
    end
        
    % run EM to approximate M
    iters = 1; MAXITERS = 2;
    OLDC = C; % old centroids matrix
    EPS = 1e-2;
    
    while iters <= MAXITERS
        
%         fprintf('on iter %d of EM...\n', iters);

%         STEP 1 - EXPECTATION
%         
%         Compute the probabilities P(x_i/b), i.e., the prob that
%         point x_i belongs to cluster core b. Use a multivariate
%         Gaussian mixture model.
        
        pjs = zeros(k,1);
        
        for i = 1:n
            for j = 1:k
                ai = A(j,1); aj = A(j,2);
                x = X(i, [ai aj]);
                mu = C(j,:);
                if sum(M(:,j) ~= 0) < 10
                    % too few points for this cluster.. 
                    % cov will bork, i.e., not positive definite.
                    pjs(j) = 0;
                else
                    pjs(j) = mvnpdf(x, mu, S(j)) * alphas(j);
                end
            end
            M(i,:) = pjs ./ sum(pjs); % memberships for the i-th point
        end
        
        % STEP 2 - MAXIMIZATION
        %
        % In this step we maximize the likelihood of the parameters
        % of the cluster cores, namely the means and covariances.
        
        alphas = sum(M, 1) ./ n; % update alpha values
        
        for i = 1:k
            ai = A(i,1); aj = A(i,2);
            Xsub = X(:, [ai aj]);
            sc = diag(M(:,i)); % scaling matrix
            XsubS = sc * Xsub;           
            C(i,:) = sum(XsubS,1) ./ sum(M(:,i)); % new centroid
           
            xc = bsxfun(@times, bsxfun(@minus, Xsub, C(i,:)), M(:,i)); % centered and scaled
            S(i) = (xc' * xc) ./ sum(M(:,i)); % new covariance matrix
        end
        
        % stop if there is no change in the centroids
        if all(abs(C(:) - OLDC(:)) < EPS) 
            break;
        end

        OLDC = C;        
        iters = iters + 1;
    end
    
    assigns = zeros(n,1);
    
    % assign each point to the cluster to which it has highest
    % pertinence.
    for i = 1:n
        [~, idx] = max(M(i,:));
        assigns(i) = idx;
    end
    
    realassigns(~noise) = assigns;
    assigns = realassigns;
    
    % now mark as noise all points that belong to clusters with less
    % than minsuprt.
    for i = 1:size(A,1)
        spt = assigns == i;
        if sum(spt) < minsuprt
            assigns(spt) = 0;
        end
    end
    
end

% find p-signatures with p = 2
% X - data set 
% S - map from attributes to sections
% P - p-signatures matrix
% A - p-signatures attributes matrix
function [P, A] = findPsignatures(X, S, minsuprt)
    d = size(X, 2); %dimensionality
    P = zeros(d * d * 10, 4); % [lower_1 upper_1 lower_2 upper_2 ; ...]
    A = zeros(d * d * 10, 2); % [att_1 att_2 ; ...]
    nused = 0; % number of used combinations
    threshold = 1e-10; % threshold for poisson acceptance
    
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
                    if ObsS > minsuprt && ObsS > ESupp && prob < threshold
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
%     [hc, edges] = histcounts(y, 'BinMethod', 'auto');
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
