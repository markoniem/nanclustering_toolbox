function D = nanmatrixdist2(X, sx, Ex, Ex2, Ex3, Ex4, idx)
% Description: 
% Compute distances between input data matrix and selected data vector. 
% Uses expected distance estimation for treating missing values. 
%
% Function call:
%          D = nanmatrixdist2(X, sx, Ex, Ex2, Ex3, Ex4, idx)
%
% Inputs:
%          X - Input data matrix
%         sx - Variance of X 
%         Ex - Expected value of X
%        Ex2 - Expected value of X.^2
%        Ex3 - Expected value of X.^3
%        Ex4 - Expected value of X.^4
%        idx - index of selected vector
%
% Output:
%         D - Calculated distances
%
Y = X(idx,:);
sy = sx(idx,:);
% Calculate omega
s = sum(sx,2) + sum(sy,2);
omega = bsxfun(@minus,X,Y);
omega = sum(omega.^2,2);
omega = omega + s;
% Calculate variance
Ey = Ex(idx,:);
Ey2 = Ex2(idx,:);
Ey3 = Ex3(idx,:);
Ey4 = Ex4(idx,:);
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
       
end

