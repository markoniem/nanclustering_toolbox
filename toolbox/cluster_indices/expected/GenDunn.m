function GDunn = GenDunn(~, Xi, sx, C, L, distance)
% Description: 
% Compute Generalized Dunn index. Algorithm uses expected distance 
% estimation for treating missing values.
%
% Function call:
%     GDunn = GenDunn(~, Xi, sx, C, L, distance)
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
%       GDunn - Value of Generalized Dunn index  
%
clusts = unique(L);
num = length(clusts);
if num == 1, GDunn = Inf; return; end
Intra = zeros(num,1);
for i = 1:num
    members = (L == clusts(i));
    Intra(i) = nancentroid([],[],nanmatrixdist(Xi(members,:), ...
        repmat(C(i,:),sum(members),1),distance,sx(members,:)),'sqe');
end
Intra = max(Intra);
Inter = min(nanpdistfunc(C,distance));
GDunn = Inter / Intra;
% Inverse
GDunn = 1 / GDunn;
              
end

