function D = ESD(X)
% Description: 
% Compute pairwise distance between observations with missing values by 
% using expected squared Euclidean distance strategy.
%
% Function call:
%        D = ESD(X)
%
% Input:
%        X - Input data set
%
% Output:
%        D - Calculated distances
%
[X, sx] = ecmnmlefunc(X);
sx = sum(sx,2);
[N, ~] = size(X);
ss = zeros(N*N-N*(N+1)/2,1);
idx = 1;
for i = 1:length(sx)-1
    for j = i+1:length(sx)
        ss(idx) = sx(i) + sx(j);
        idx = idx + 1;
    end
end
D = sqrt((pdist(X)').^2+ss);

end

