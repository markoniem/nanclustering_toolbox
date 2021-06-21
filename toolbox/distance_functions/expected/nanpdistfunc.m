function D = nanpdistfunc(X, distance, sx)
% Description: 
% Compute pairwise distances between pairs of observations in data 
% matrix. Uses expected distance estimation for treating missing values. 
%
% Function call:
%        D = nanpdistfunc(X, distance, sx)
%
% Input:
%         X - Input data matrix 
%  distance - Selected distance metric 
%             Alternatives: 
%             'euc' - Euclidean distance 
%             'sqe' - squared Euclidean distance
%       sx - Variance of matrix 
%
% Output:
%        D - Distances between observations
%
if nargin == 2
    sx = zeros(size(X));
end
switch distance
    case 'euc'
        % Calculate omega
        [N, ~] = size(X);
        s = sum(sx,2);
        ss = zeros(N*N-N*(N+1)/2,1);
        idx = 1;
        for i = 1:length(s)-1
            for j = i+1:length(s)
                ss(idx) = s(i) + s(j);
                idx = idx + 1;
            end
        end
        omega = (pdist(X)').^2 + ss;
        % Calculate variance
        Ex = X;
        Ex2 = X.^2 + sx;
        Ex3 = X.^3 + 3*X.*sx;
        Ex4 = X.^4 + 6*(X.^2).*sx + 3*sx.^2;
        var1 = zeros(N*N-N*(N+1)/2,1);
        var2 = zeros(N*N-N*(N+1)/2,1);
        idx = 1;
        for i = 1:length(s)-1
            for j = i+1:length(s)
                var1(idx) = sum(Ex4(i,:) + Ex4(j,:) - 4*Ex3(i,:).*Ex(j,:) - ... 
                    4*Ex(i,:).*Ex3(j,:) + 6*Ex2(i,:).*Ex2(j,:),2);
                var2(idx) = sum((Ex2(i,:) - 2*Ex(i,:).*Ex(j,:) + Ex2(j,:)).^2,2);
                idx = idx + 1;
            end
        end
        var = var1 - var2;
        var(var<0.0000001) = 0;
        % Calculate EED
        m = (omega.^2)./var;
        D = exp(gammaln(m+0.5) - gammaln(m));
        D = D.*((omega./m).^(0.5));
        ind = isnan(D);
        D(ind) = sqrt(omega(ind));   
    case 'sqe'
        % Calculate ESD
        [N, ~] = size(X);
        s = sum(sx,2);
        ss = zeros(N*N-N*(N+1)/2,1);
        idx = 1;
        for i = 1:length(s)-1
            for j = i+1:length(s)
                ss(idx) = s(i) + s(j);
                idx = idx + 1;
            end
        end
        D = (pdist(X)').^2 + ss;    
end

end

