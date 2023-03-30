function error = rmsefunc(X1, X2, M, N)
% Description: 
% Computes root mean square error between real and estimated values of
% distances.
%
% Inputs:
%
%         X1 - Real distances between complete observations
%         X2 - Estimated distances between observations with missing values
%          M - Number of incomplete observations
%          N - Number of all observation
%
% Output:
%
%      error - Root mean square error
%
alpha = M*N - M*(M + 1) / 2;
error = sqrt((1/alpha)*sum((X1 - X2).^2));
end

