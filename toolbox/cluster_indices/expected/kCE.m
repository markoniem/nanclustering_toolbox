function kCEVal = kCE(~, Xi, sx, C, L, ~)
% Description: 
% Compute kCE-index. Algorithm uses expected distance 
% estimation for treating missing values. 
%
% Function call:
%    kCEVal = kCE(~, Xi, sx, C, L, ~)
%
% Inputs:
%        Xi - Imputed data set
%        sx - Variance of missing values in original data set 
%         C - Matrix of cluster centroids 
%         L - Cluster labels for each observation
%
% Output:
%       kCEVal - Value of kCE-index 
%
k = size(C,1);
kCEVal = k*nansum(nandistfuncp2(C,Xi,L,'sqe',sx),1);

end

