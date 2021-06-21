function Silh = Silhouette(~, Xi, sx, ~, L, ~)
% Description: 
% Compute Silhouette index. Algorithm uses expected distance 
% estimation for treating missing values.
%
% Function call:
%      Silh = Silhouette(~, Xi, sx, ~, L, ~)
%
% Inputs:
%        Xi - Imputed data set
%        sx - Variance of missing values in original data set  
%         L - Cluster labels for each observation
%
% Output:
%       Silh - Value of Silhouette index    
%
cnames = unique(L);
k = length(cnames);
if k == 1, Silh = Inf; return; end
n = length(L);
mbrs = (repmat(1:k,n,1) == repmat(L,1,k));
avgDWithin = Inf(n,1);
avgDBetween = Inf(n,k);
Ex = Xi;
Ex2 = Xi.^2 + sx;
Ex3 = Xi.^3 + 3*Xi.*sx;
Ex4 = Xi.^4 + 6*(Xi.^2).*sx + 3*sx.^2;
for j = 1:n
    distj = nansum(nanmatrixdist2(Xi,sx,Ex,Ex2,Ex3,Ex4,j),2);
    for i = 1:k
        if i == L(j)
            mbrs1 = mbrs;
            mbrs1(j,i) = 0;
            avgDWithin(j) = nancentroid([],[],distj(mbrs1(:,i)),'sqe');
        else
            avgDBetween(j,i) = nancentroid([],[],distj(mbrs(:,i)),'sqe');
        end
    end
end
minavgDBetween = min(avgDBetween, [], 2);
Silh = (minavgDBetween - avgDWithin) ./ max(avgDWithin,minavgDBetween);
Silh = nancentroid([],[],Silh,'sqe');
% Inverse
Silh = 1 / Silh;

end

