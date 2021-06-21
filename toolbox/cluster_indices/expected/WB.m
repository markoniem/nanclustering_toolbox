function WBVal = WB(Xm, Xi, sx, C, L, distance)
% Description: 
% Compute WB-index. Algorithm uses expected distance estimation for 
% treating missing values. 
%
% Function call:
%     WBVal = WB(Xm, Xi, sx, C, L, distance)
%
% Inputs:
%        Xm - Original data set with missing values
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
%       WBVal - Value of WB-index
%
K = size(C,1);
if K == 1, WBVal = Inf; return; end
Ni = zeros(K,1);
for i = 1:K, l = find(L==i); Ni(i) = length(l); end
Intra = nansum(nandistfuncp2(C,Xi,L,'sqe',sx));
Inter = nansum(Ni.*nansum(nanmatrixdist(C,nancentroid(Xi,sx,Xm,distance),'sqe'),2));
WBVal = K*Intra / Inter;

end

