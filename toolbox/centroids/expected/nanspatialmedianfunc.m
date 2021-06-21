function C = nanspatialmedianfunc(X, sx, X_org, varargin)  
% Description: 
% Estimate spatial median value of input data set. Algorithm uses expected 
% distance estimation for treating missing values. 
%
% Function calls:
%         C = nanspatialmedianfunc(X, sx, X_org)
%         C = nanspatialmedianfunc(X, sx, X_org, max_iter)
%         C = nanspatialmedianfunc(X, sx, X_org, max_iter, tol)
%
% Inputs: 
%         X - Imputed data set
%        sx - Variance of missing values in original data set
%     X_org - Original data set
%  max_iter - Maximum number of iterations. Default value: 100 
%       tol - Algorithm stopping criterion. Default value: 1e-5
%
% Output: 
%         C - Spatial median value of data set
%
defaultMaxIter = 100;
defaultTol = 1e-5; 
%
p = inputParser;
addOptional(p, 'max_iter', defaultMaxIter, @(x) isnumeric(x));
addOptional(p, 'tol', defaultTol, @(x) isnumeric(x));
parse(p,varargin{:});
max_iter = p.Results.max_iter;
tol = p.Results.tol;
%
u = zeros(1,size(X_org,2));
for i = 1:size(X_org,2)
    I = ~isnan(X_org(:,i));
    u(i) = median(X_org(I,i));
end

P = ~isnan(X_org);
iters = 0;
w = 1.5;
while iters < max_iter
    iters = iters + 1;
    D = eucdistfunc(u,X,ones(size(X,1),1),sx);
    a = 1./(D+sqrt(sqrt(eps)));
    ax = nansum(bsxfun(@times,a,X_org));
    a = bsxfun(@times,P,a);
    v = (1./sum(a)).*ax;
    u1 = u + w*(v-u);
    if norm(u1-u,inf) < tol
        break;
    end
    u = u1;
end
C = u1;

end

function D = eucdistfunc(C, Xi, L, sx) 
% Description: Calculate euclidean distances between observation and 
% centroid matrices. Uses expected distance estimation for treating 
% missing values. 
% Input:
%        C - Matrix of cluster centroids
%       Xi - Inputed data set
%        L - Cluster labels for each observation
%       sx - Variance of missing values in original data set
% Output:
%        D - Distance matrix   
sx = sum(sx,2);
D = bsxfun(@minus,Xi,C(L,:));
D = sqrt(sum(D.^2,2) + sx);

end

