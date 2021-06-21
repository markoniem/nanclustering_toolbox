function CH = CalinskiHarabasz(Xm, Xi, sx, C, L, distance)
% Description: 
% Compute Calinski-Harabasz index. Algorithm uses expected distance 
% estimation for treating missing values.
%
% Function call:
%        CH = CalinskiHarabasz(Xm, Xi, sx, C, L, distance)
%
% Inputs:
%        Xm - Input data set consisting missing values 
%        Xi - Imputed data set
%        sx - Variance of missing values in Xm 
%         C - Matrix of cluster centroids 
%         L - Cluster labels for each observation
%  distance - Selected distance metric 
%             Alternatives: 
%             'euc' - Euclidean distance 
%             'sqe' - squared Euclidean distance
%
% Output:
%       CH - Value of Calinski-Harabasz index    
%
K = size(C,1);
if K == 1, CH = Inf; return; end
N = size(Xm,1);
Ni = zeros(K,1);
for i = 1:K, l = find(L==i); Ni(i) = length(l); end
Intra = nansum(nandistfuncp2(C,Xi,L,'sqe',sx));
Inter = nansum(Ni.*nanmatrixdist(C,nancentroid(Xi,sx,Xm,distance),'sqe'));
CH = ((N-K)/(K-1))*(Inter/Intra);
% Inverse
CH = 1 / CH;

end

