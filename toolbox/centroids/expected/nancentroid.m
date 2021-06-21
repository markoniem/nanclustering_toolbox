function C = nancentroid(Xi, sx, X, distance) 
% Description: 
% Estimate centroid value of input data set. Algorithm uses expected 
% distance estimation for treating missing values.
%
% Function call:
%         C = nancentroid(Xi, sx, X, distance)   
%
% Inputs: 
%         Xi - Imputed data set
%         sx - Variance of original data set with missing values
%          X - Original data set 
%   distance - Selected norm distance metric 
%              Alternatives: 
%              'euc' - Euclidean distance
%              'sqe' - squared Euclidean distance
%
% Output: 
%         C - Centroid value of data set 
%
switch distance
    case 'euc'
        % Spatial median value of data
        C = nanspatialmedianfunc(Xi,sx,X);
    case 'sqe'
        nan = isnan(X);
        X(nan) = 0;
        % Mean value of data
        C(1,:) = sum(X,1) ./ sum(~nan,1);  
end

end

