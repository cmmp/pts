  
%     MINPTS = 20;
%     K = 0;
%     toRemove = logical(zeros(size(A,1),1)); % remove non-significant p-sigs.
%     for i = 1:size(A,1)
%         ai = A(i,1); aj = A(i,2);
%         pi = (X(:,ai) >= P(i,1)) & (X(:,ai) <= P(i,2));
%         pj = (X(:,aj) >= P(i,3)) & (X(:,aj) <= P(i,4));
%         pmix = find(pi & pj);
%         Xsub = X(pmix, [ai aj]);
%         Eps = epsestimate(Xsub, MINPTS);
%         [clazz, ~] = dbscan(Xsub, MINPTS, Eps);
%         if all(clazz == -1)
%             toRemove(i) = 1;
%         end
%     end
%    A(toRemove,:) = [];
    
%     assigns = zeros(n,1);
    
    % the first order of business is splitting clusters that are
    % incorrectly together
%     
%     MINPTS = 10;
%     K = 0;
%     B = zeros(size(A,1) * 10, 2);    
%     
%     q = 1;
%     for i = 1:size(A,1)
%         ai = A(i,1); aj = A(i,2);
%         pi = (X(:,ai) >= P(i,1)) & (X(:,ai) <= P(i,2));
%         pj = (X(:,aj) >= P(i,3)) & (X(:,aj) <= P(i,4));
%         pmix = find(pi & pj);
%         
%         if isempty(pmix)
%             error('got an empty cluster at dbscan split phase! die.');
%         end
%         
%         % now use dbscan on pmix to separate any clusters that might
%         % have been erroneously merged
%         Xsub = X(pmix, [ai aj]);
%         Eps = epsestimate(Xsub, MINPTS);
%         [clazz, ~] = dbscan(Xsub, MINPTS, Eps);
%         clazz(clazz == -1) = 0; % noise is 0 to me
%         idx = find(clazz ~= 0);
%         clazz(idx) = clazz(idx) + K;
%         uclass = length(unique(clazz(clazz ~= 0)));
%         B(q:(q+uclass-1),:) = repmat(A(i,:), uclass, 1);
%         q = q + uclass;
%         K = max(clazz);
%         assigns(pmix) = clazz;
%     end
    
%       clf(); plotWindow(X, assigns);
    
%     A = B(B(:,1) > 0,:); % remove unused spaces    
%     B = 0; % just free some space
    
%     % we should remove clusters that have no points
%     toRemove = logical(zeros(size(A,1), 1));
%     for i = 1:size(A,1)
%         if sum(assigns == i) == 0
%             toRemove(i) = 1;
%         end
%     end
%     A(toRemove,:) = [];