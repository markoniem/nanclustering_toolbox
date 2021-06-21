function D = nanmatrixdist(X, Y, distance, sx, sy)
% Description: 
% Compute distance between two equal sized matrices. Uses expected 
% distance estimation for treating missing values. 
%
% Function call:
%          D = nanmatrixdist(X, Y, distance, sx, sy)
%
% Input:
%          X - First data set
%          Y - Second data set  
%   distance - Selected distance metric 
%              Alternatives: 
%              'euc' - Euclidean distance 
%              'sqe' - squared Euclidean distance
%         sx - Variance of X
%         sy - Variance of Y
%
% Output:
%         D - Distance matrix 
%
if nargin == 3
    sx = zeros(size(X));
    sy = zeros(size(Y));
end
if nargin == 4
    sy = zeros(size(Y));
end
switch distance
    case 'euc'    
        % Calculate omega
        s = sum(sx,2) + sum(sy,2);
        omega = bsxfun(@minus,X,Y);       
        omega = sum(omega.^2,2);
        omega = omega + s;
        % Calculate variance
        Ex = X;
        Ex2 = X.^2 + sx;
        Ex3 = X.^3 + 3*X.*sx;
        Ex4 = X.^4 + 6*(X.^2).*sx + 3*sx.^2;
        Ey = Y;
        Ey2 = Y.^2 + sy;
        Ey3 = Y.^3 + 3*Y.*sy;
        Ey4 = Y.^4 + 6*(Y.^2).*sy + 3*sy.^2;        
        var = sum(bsxfun(@plus,Ex4,Ey4) - 4*bsxfun(@times,Ex3,Ey) - ...
                    4*bsxfun(@times,Ex,Ey3) + 6*bsxfun(@times,Ex2,Ey2),2) - ...
                        sum((bsxfun(@plus,Ex2-2*bsxfun(@times,Ex,Ey),Ey2)).^2,2);
        var(var<0.0000001) = 0;
        % Calculate EED
        m = (omega.^2)./var;
        D = exp(gammaln(m + 0.5) - gammaln(m));
        D = D.*((omega./m).^(0.5));
        ind = isnan(D);
        D(ind) = sqrt(omega(ind));
    case 'sqe'
        % Calculate ESD
        sx = sum(sx,2); 
        sy = sum(sy,2);
        s = sx + sy;  
        D = bsxfun(@minus,X,Y);       
        D = sum(D.^2,2);
        D = D + s;            
end
             
end

