function L = nandistfuncp1(C, X, distance, sx) 
% Description: 
% Compute nearest centroid for each observation. Uses expected distance 
% estimation for treating missing values. 
%
% Function call:
%         L = nandistfuncp1(C, X, distance, sx)
%
% Inputs:  
%         C - Matrix of cluster centroids
%         X - Input data set
%  distance - Selected distance metric 
%             Alternatives: 
%             'euc' - Euclidean distance 
%             'sqe' - squared Euclidean distance
%        sx - Variance of input data set
%
% Output: 
%         L - Cluster labels for each observation
%
if nargin == 3
    sx = zeros(size(X));
end
switch distance
    case 'euc'
        % Calculate EED
        D = pdist2(X,C).^2;
        sx = sum(sx,2);
        D = bsxfun(@plus,D,sx);
        [~, L] = min(sqrt(D),[],2);
        L = L(:);
    case 'sqe'
        % Calculate ESD
        D = pdist2(X,C).^2;
        sx = sum(sx,2);
        D = bsxfun(@plus,D,sx);
        [~, L] = min(D,[],2);
        L = L(:);
end

end

