function [D, success] = IST_MC(X)
% Description: 
% Compute pairwise distance between observations with missing values 
% by using matrix completion via iterative soft-thresholding.
%
% Function call:
%        D = IST_MC(X)
%
% Input:
%        X - Input data set
%
% Output:
%        D - Calculated distances
%  success - Indicator flag for successfully finished algorithm  
%
success = 1;
sizeX = size(X); 
x = X(:);
IDX = ~isnan(x);
M = zeros(length(x),1);
M(IDX) = 1;
x(~IDX) = 0;
%
err = 1e-6;
x_initial = zeros(prod(sizeX),1);
insweep = 200;
tol = 1e-4;    
decfac = 0.9;
%
y = M.*x;
x = x_initial;
lambdaInit = decfac*max(abs(M.*y)); lambda = lambdaInit;
f_current = norm(y-M.*x) + lambda*norm(x,1);

while lambda > lambdaInit*tol
    for ins = 1:insweep
        f_previous = f_current;
        x = x + (y - M.*x);
        [U,S,V] = svd(reshape(x,sizeX),'econ');
        s = sign(diag(S)).*max(0,abs(diag(S))-lambda/2);
        S = diag(s);
        X = U*S*V';
        x = X(:);
        f_current = norm(y-M.*x) + lambda*norm(x,1);
        if norm(f_current-f_previous)/norm(f_current + f_previous)<tol
            break;
        end
    end
    if norm(y-M.*x)<err
        break;
    end
    lambda = decfac*lambda;
end
D = pdist(X)';

end

