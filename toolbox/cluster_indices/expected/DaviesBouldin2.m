function DB = DaviesBouldin2(~, Xi, sx, C, L, distance)
% Description: 
% Compute Davies-Bouldin* index. Algorithm uses expected distance 
% estimation for treating missing values.
%
% Function call:
%        DB = DaviesBouldin2(~, Xi, sx, C, L, distance)
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
%       DB - Value of Davies-Bouldin* index    
clusts = unique(L);
num = length(clusts);
if num == 1, DB = Inf; return; end
aveWithinD = zeros(num,1);
for i = 1:num
    members = (L == clusts(i));
    aveWithinD(i) = nancentroid([],[],nanmatrixdist(Xi(members,:), ...
        repmat(C(i,:),sum(members),1),distance,sx(members,:)),'sqe');
end
interD = nandistfunc(C,C,distance);
interD(logical(eye(size(interD)))) = Inf;
R = zeros(num);
cnt = 1;
for i = 1:num
    for j = i+1:num
        R(i,j) = (aveWithinD(i)+aveWithinD(j));
        cnt = cnt + 1;
    end
end
R = R + R';
RI = max(R,[],1) ./ min(interD,[],1);
DB = nancentroid([],[],RI(:),'sqe');
       
end

