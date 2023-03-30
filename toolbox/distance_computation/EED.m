function [D, success] = EED(X)
% Description: 
% Computes pairwise distance between observations with missing values 
% using expected Euclidean distance strategy.
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
D = D(B==1);

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
%
% Compute omega
omega = pdist2(X,Y).^2;
omega = bsxfun(@plus,omega,sum(sx,2));
omega = bsxfun(@plus,omega,sum(sy,2)');
% Compute variance
Ex = X;
Ex2 = X.^2 + sx;
Ex3 = X.^3 + 3*X.*sx;
Ex4 = X.^4 + 6*(X.^2).*sx + 3*sx.^2;
Ey = Y;
Ey2 = Y.^2 + sy;
Ey3 = Y.^3 + 3*Y.*sy;
Ey4 = Y.^4 + 6*(Y.^2).*sy + 3*sy.^2;
I = ones(size(Y));
I2 = ones(size(X));
var = Ex4*I' + I2*Ey4' - 4*Ex3*Ey' - 4*Ex*Ey3' + 6*Ex2*Ey2' - ...
    (Ex2*I' - 2*Ex*Ey' + I2*Ey2').^2;
var(var<0.0000001) = 0;
% Compute EED
m = (omega.^2)./var;
D = exp(gammaln(m+0.5) - gammaln(m));
D = D.*((omega./m).^(0.5));
ind = isnan(D);
D(ind) = sqrt(omega(ind));

end


