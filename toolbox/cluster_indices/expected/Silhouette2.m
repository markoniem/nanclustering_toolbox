function Silh = Silhouette2(~, Xi, sx, C, L, distance)
% Description: 
% Compute Silhouette* index. Algorithm uses expected distance estimation for 
% treating missing values and it is speeded-up version of original Silhouette 
% index. More specify, distances are computed between observations 
% and nearest centroids instead between observations and another observations.
%
% Function call:
%      Silh = Silhouette2(~, Xi, sx, C, L, distance)
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
%       Silh - Value of Silhouette* index   
%
clusts = unique(L);
num = length(clusts);
if num == 1, Silh = Inf; return; end
avgDWithin = Inf(num,1);
avgDBetween = Inf(num,num);
for i = 1:num
    members = (L == clusts(i));
    for j = 1:num
        if j==i
            avgDWithin(i) = nancentroid([],[],nanmatrixdist(Xi(members,:), ...
                repmat(C(j,:),sum(members),1),distance,sx(members,:)),'sqe');
        else
            avgDBetween(i,j) = nancentroid([],[],nanmatrixdist(Xi(members,:), ...
                repmat(C(j,:),sum(members),1),distance,sx(members,:)),'sqe');
        end  
    end  
end
minavgDBetween = min(avgDBetween, [], 2);
Silh = (minavgDBetween - avgDWithin) ./ minavgDBetween;
Silh = nancentroid([],[],Silh,'sqe');
% Inverse
Silh = 1 / Silh;

end

