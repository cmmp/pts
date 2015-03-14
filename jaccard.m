% Compute the jaccard coefficient in [0,1], 
% bigger values are better.
% 
% Usage:
% j = jaccard(labels, expect)
%
% Arguments:
% labels - array of cluster labels
% expect - supervision array of the same length
% 
% Returns:
% j - jaccard coefficient in range [0,1]
function j = jaccard(labels, expect)
    a = 0; b = 0; c = 0;
    j = 0;
    if ~isequal(size(labels), size(expect))
        error('labels and expect must have the same length!')
    end
    
    N = length(labels);
    
    for i = 1:(N-1)
        for j = (i+1):N
            if labels(i) == labels(j) && expect(i) == expect(j)
                a = a + 1;
            elseif labels(i) ~= labels(j) && expect(i) == expect(j)
                b = b + 1;
            elseif labels(i) == labels(j) && expect(i) ~= expect(j)
                c = c + 1;
            end
        end
    end
    
    j = a / (a + b + c);
    
end