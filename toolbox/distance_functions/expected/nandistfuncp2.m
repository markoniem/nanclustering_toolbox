function D = nandistfuncp2(C, X, L, distance, sx) 
% Description: 
% Compute distances between two matrices. Uses expected distance 
% estimation for treating missing values. 
%
% Function call:
%        D = nandistfuncp2(C, X, L, distance, sx)
%
% Inputs:
%        C - Matrix of cluster centroids 
%        X - Input data set 
%        L - Cluster labels for each observation
% distance - Selected distance metric 
%            Alternatives: 
%            'euc' - Euclidean distance 
%            'sqe' - squared Euclidean distance
%       sx - Variance of input data set
%
% Output:
%        D - Distances to nearest centroids
%
if nargin == 4
    sx = zeros(size(X));
end
switch distance
    case 'euc'
        % Calculate omega
        omega = bsxfun(@minus,X,C(L,:));
        omega = sum(omega.^2,2) + sum(sx,2);  
        % Calculate variance
        Ex = X;
        Ex2 = X.^2 + sx;
        Ex3 = X.^3 + 3*X.*sx;
        Ex4 = X.^4 + 6*(X.^2).*sx + 3*sx.^2;
        Y = C(L,:);
        Ey = Y;
        Ey2 = Y.^2;
        Ey3 = Y.^3;
        Ey4 = Y.^4;
        var = sum(Ex4 + Ey4 - 4*Ex3.*Ey - 4*Ex.*Ey3 + 6*Ex2.*Ey2,2) - ...
            sum((Ex2 - 2*Ex.*Ey + Ey2).^2,2);
        var(var<0.0000001) = 0; 
        % Calculate EED
        m = (omega.^2)./var;
        D = exp(gammaln(m+0.5)-gammaln(m));
        D = D.*((omega./m).^(0.5));
        ind = isnan(D);
        D(ind) = sqrt(omega(ind));
    case 'sqe'
        % Calculate ESD
        sx = sum(sx,2);
        D = bsxfun(@minus,X,C(L,:));
        D = sum(D.^2,2) + sx;
end
        
end

