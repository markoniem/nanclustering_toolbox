function PBMVal = PBM(Xm, Xi, sx, C, L, distance)
% Description: 
% Compute Pakhira-Bandyopadhyay-Maulik index. Algorithm uses expected 
% distance estimation for treating missing values. 
%
% Function call:
%    PBMVal = PBM(Xm, Xi, sx, C, L, distance)
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
%       PBMVal - Value of Pakhira-Bandyopadhyay-Maulik index
%
k = size(C,1);
if k == 1, PBMVal = Inf; return; end
Intra = nansum(nanmatrixdist(Xi,C(L,:),distance,sx));
Inter = max(nanpdistfunc(C,distance));
err = nansum(nanmatrixdist(Xi,nancentroid(Xi,sx,Xm,distance),distance,sx));
PBMVal = (Inter*err)/(k*Intra);
% Inverse
PBMVal = 1 / PBMVal;
 
end

