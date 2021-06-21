function D = nandistfuncp2(C, X, L, distance) 
% Description: 
% Estimate distances between observation and centroid matrices. Uses 
% partial distance strategy for treating missing values.
%
% Function call:
%         D = nandistfuncp2(C, X, L, distance)
%
% Input:
%         C - Matrix of cluster centroids
%         X - Input data set
%         L - Cluster labels for each observation
%  distance - Selected distance metric 
%             Alternatives: 
%             'sqe' - squared Euclidean distance 
%             'euc' - Euclidean distance 
%             'cit' - City block distance  
%            'scit' - squared City block distance  
%
% Output:
%         D -  Distances to nearest centroids
%
switch distance
    case 'sqe'
        D = bsxfun(@minus,X,C(L,:));
        D = nansum(D.^2,2);
        D = D.*(size(X,2)./sum(~isnan(X),2));
    case 'euc'
        D = bsxfun(@minus,X,C(L,:));
        D = nansum(D.^2,2);
        D = D.*(size(X,2)./sum(~isnan(X),2));
        D = sqrt(D);
    case 'cit'
        D = bsxfun(@minus,X,C(L,:));
        D = nansum(abs(D),2);
        D = D.*(size(X,2)./sum(~isnan(X),2));
    case 'sqcit'
        D = bsxfun(@minus,X,C(L,:));
        D = nansum(abs(D),2);
        D = D.*(size(X,2)./sum(~isnan(X),2));
        D = D.^2;
end
        
end

