function RT = RayTuri(~, Xi, sx, C, L, distance)
% Description: 
% Compute Ray-Turi index. Algorithm uses expected distance 
% estimation for treating missing values.
%
% Function call:
%        RT = RayTuri(~, Xi, sx, C, L, distance)
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
%       RT - Value of Ray-Turi index
%
K = size(C,1);
if K == 1, RT = Inf; return; end
Intra = nansum(nandistfuncp2(C,Xi,L,distance,sx));
Inter = min(nanpdistfunc(C,distance));
RT = Intra / Inter;

end

