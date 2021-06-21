function WG = WemmertGancarski(~, Xi, sx, C, L, distance)
% Description: 
% Compute Wemmert-Gancarski index. Algorithm uses expected distance 
% estimation for treating missing values.
%
% Function call:
%        WG = WemmertGancarski(~, Xi, sx, C, L, distance)
%
% Inputs:
%        Xi - Imputed data set
%        sx - Variance of missing values in original data set 
%         C - Matrix of cluster centroids 
%         L - Cluster labels for each observation
%  distance - Selected distance metric 
%             Alternatives: 
%             'euc' - Euclidean distance 
%             'sqe' - squared Euclidean distance
%
% Output:
%       WG - Value of Wemmert-Gancarski index
%
k = size(C,1);
if k == 1, WG = Inf; return; end
n = size(C,2);
Inter = 0;
C2 = C;
for i=1:k
    C = C2;
    I = find(L == i);
    C(i,:) = realmax/(10^6)*ones(1,n);
    OthClustdists = min(nandistfunc(Xi(I,:),C,distance,sx(I,:)),[],2);
    RM = nanmatrixdist(Xi(I,:),repmat(C2(i,:),length(I),1),distance, ...
                            sx(I,:))./OthClustdists;
    Inter = Inter + length(I) - nansum(RM);
end
WG = Inter;
% Inverse
WG = 1 / WG;

end

