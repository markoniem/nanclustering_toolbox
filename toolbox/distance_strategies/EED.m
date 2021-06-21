function D = EED(X)
% Description: 
% Compute pairwise distance between observations with missing values 
% by using expected Euclidean distance strategy.
%
% Function call:
%        D = EED(X)
%
% Input:
%        X - Input data set
%
% Output:
%        D - Calculated distances
%
[X, sx] = ecmnmlefunc(X);
[N, ~] = size(X);
ss = zeros(N*N-N*(N+1)/2,1);
idx = 1;
for i = 1:size(sx,1)-1
    for j = i+1:size(sx,1)
        ss(idx) = sum(sx(i,:),2) + sum(sx(j,:),2);
        idx = idx + 1;
    end
end
omega = (pdist(X)').^2+ss;
Ex = X;
Ex2 = X.^2 + sx;
Ex3 = X.^3 + 3*X.*sx;
Ex4 = X.^4 + 6*(X.^2).*sx + 3*sx.^2;
var1 = zeros(N*N-N*(N+1)/2,1);
var2 = zeros(N*N-N*(N+1)/2,1);
idx = 1;
for i = 1:size(X,1)-1
    for j = i+1:size(X,1)
        var1(idx) = sum(Ex4(i,:) + Ex4(j,:) - 4*Ex3(i,:).*Ex(j,:) - ... 
                        4*Ex(i,:).*Ex3(j,:) + 6*Ex2(i,:).*Ex2(j,:),2);  
        var2(idx) = sum((Ex2(i,:) - 2*Ex(i,:).*Ex(j,:) + ... 
                        Ex2(j,:)).^2,2);
        idx = idx + 1;
    end
end
var = var1 - var2;
var(var<0.0000001) = 0;
m = (omega.^2)./var; 
D = exp(gammaln(m+0.5) - gammaln(m));
D = D.*((omega./m).^(0.5));
ind = isnan(D);
D(ind) = sqrt(omega(ind));

end

