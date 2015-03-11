
% p3c function
function p3c(X, sup)
    % identify relevant sections in attributes
    S = findSections(X);
end

function S = findSections(X)
    [N, d] = size(X);
    S = containers.Map('KeyType', 'int64', 'ValueType', 'any');
    
    for i = 1:d
        S(i) = analyzeAttribute(X(:,i));
    end

    %disp('cluster cores phase')
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
