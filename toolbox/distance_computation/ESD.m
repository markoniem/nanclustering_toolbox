function [D, success] = ESD(X)
% Description: 
% Computes pairwise distance between observations with missing values 
% using expected squared Euclidean distance strategy.
%
% Input:
%        X - Input data set with missing values
%
% Output:
%        D - Computed distances
%  success - Indicator flag for successfully finished algorithm  
%
[Xi, sx, success] = ecmnmlefunc(X);
D = nandistfunc(Xi, Xi, sx, sx);
N = size(X,1);
B = tril(ones(N,N),-1);
D = sqrt(D(B==1));

end

function D = nandistfunc(X, Y, sx, sy)
% Description: 
% Computes distances between two matrices. Uses expected distance 
% estimation for treating missing values. 
%
% Inputs: 
%         X - Imputed data set using conditional mean
%         Y - Second imputed data set using conditional mean
%
% Output:
%       D - Distance matrix
%
if nargin == 3
    sx = zeros(size(X));
    sy = zeros(size(Y));
end
if nargin == 4
    sy = zeros(size(Y));
end
% Compute ESD
D = pdist2(X,Y).^2;
sx = sum(sx,2);
sy = sum(sy,2);
D = bsxfun(@plus,D,sx);
D = bsxfun(@plus,D,sy');

        
end

