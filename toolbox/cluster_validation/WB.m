function WBVal = WB(X, C, L, weights, Xi, sx, varargin)
% Description: 
% Computes WB-index value. Three options are available 
% for treating missing values. See documentation. 
%
% Inputs:
%               X - Input data set with missing values
%               C - Cluster centroids 
%               L - Cluster labels 
%         weights - Weights of observations. Default is a vector of ones.
%              Xi - Imputed data using conditional mean imputation
%              sx - Variance of data set
%
%  'TreatMissing' - Missing values treating strategy. Default is 'ads'.      
%                   Alternatives:
%                   'ads' - available data strategy
%                   'pds' - partial distance strategy
%                   'exp' - expected distance strategy                      
%      'Distance' - Selected distance metric. Default is 'euc'.
%                   Alternatives: 
%                   'sqe' - squared Euclidean distance
%                   'euc' - Euclidean distance
%                   'cit' - City block distance 
%                   Note: Expected City block distance is not supported.
%
% Output:
%        WB - Value of WB-index 
%
pnames = {'treatmissing' 'distance'};
dflts =  {'ads' 'euc'};
[treatmissing, distance] = internal.stats.parseArgs(pnames, dflts, varargin{:});
treatmissingNames = {'ads','pds','exp'};
distNames = {'sqe','euc','cit'};
if isnumeric(treatmissing) || ~ismember(treatmissing,treatmissingNames)
    error('Invalid treatmissing');
end
if isnumeric(distance) || ~ismember(distance,distNames)
    error('Invalid distance');
end
if (strcmp(treatmissing,'exp') && strcmp(distance,'cit'))
    error('Expected City block is not supported');
end
%
if strcmp(distance,'sqe')
    nancentroid = @nanmean;
elseif strcmp(distance,'cit')
    nancentroid = @nanmedian;
end
if ~strcmp(treatmissing,'exp')
    dist2 = distance;
    if strcmp(distance,'euc'), dist2 = 'sqe'; end
    if strcmp(distance,'cit'), dist2 = 'sqcit'; end
    if strcmp(treatmissing,'ads')
        nandistfuncp2 = @nandistfuncp2_ads;
        nanmatrixdist = @nanmatrixdist_ads;
        if strcmp(distance,'euc')
            nancentroid = @nanspatialmedian_ads;
        end
    elseif strcmp(treatmissing,'pds')
        nandistfuncp2 = @nandistfuncp2_pds;
        nanmatrixdist = @nanmatrixdist_pds;
        if strcmp(distance,'euc')
            nancentroid = @nanspatialmedian_pds;
        end      
    end    
    K = size(C,1);
    if K == 1, WBVal = Inf; return; end
    Ni = zeros(K,1);
    for i = 1:K, l = (L==i); Ni(i) = sum(weights(l)); end
    Intra = nansum(weights.*nandistfuncp2(C,X,L,dist2));
    Inter = nansum(Ni.*nansum(nanmatrixdist(C,nancentroid(X,weights),dist2),2));
    WBVal = K*Intra / Inter;
else   
    if strcmp(distance,'euc')
        nancentroid = @nanspatialmedian_exp;
    end    
    K = size(C,1);
    if K == 1, WBVal = Inf; return; end
    Ni = zeros(K,1);
    for i = 1:K, l = (L==i); Ni(i) = sum(weights(l)); end
    Intra = nansum(weights.*nandistfuncp2_exp(C,Xi,L,'sqe',sx));
    Inter = nansum(Ni.*nansum(nanmatrixdist_exp(C,nancentroid(X,weights,Xi,sx),'sqe'),2));
    WBVal = K*Intra / Inter;
end

end

function D = nandistfuncp2_ads(C, X, L, distance) 
% Description: 
% Computes distances between observations with missing values and nearest 
% centroids. Uses available data strategy for treating missing values. 
%
% Inputs:
%         C - Cluster centroids
%         X - Input data set
%         L - Cluster labels 
%  distance - Selected distance metric 
%             Alternatives: 
%             'sqe' - squared Euclidean distance
%             'euc' - Euclidean distance 
%             'cit' - City block distance 
%             'sqcit' - squared City block distance 
%
% Output:
%         D -  Distances to nearest centroids 
%
switch distance
    case 'sqe'
        D = bsxfun(@minus,X,C(L,:));
        D = nansum(D.^2,2);
    case 'euc'
        D = bsxfun(@minus,X,C(L,:));
        D = nansum(D.^2,2);
        D = sqrt(D);        
    case 'cit'
        D = bsxfun(@minus,X,C(L,:));
        D = nansum(abs(D),2);
    case 'sqcit'
        D = bsxfun(@minus,X,C(L,:));
        D = nansum(abs(D).^2,2);        
end
        
end

function D = nandistfuncp2_pds(C, X, L, distance) 
% Description: 
% Computes distances between observations with missing values and nearest 
% centroids. Uses partial distance strategy for treating missing values. 
%
% Input:
%         C - Cluster centroids
%         X - Input data set
%         L - Cluster labels
%  distance - Selected distance metric 
%             Alternatives: 
%             'sqe' - squared Euclidean distance 
%             'euc' - Euclidean distance 
%             'cit' - City block distance  
%            'scit' - squared City block distance  
%
% Output:
%         D -  Distances to nearest centroids
%
switch distance
    case 'sqe'
        D = bsxfun(@minus,X,C(L,:));
        D = nansum(D.^2,2);
        D = D.*(size(X,2)./sum(~isnan(X),2));
    case 'euc'
        D = bsxfun(@minus,X,C(L,:));
        D = nansum(D.^2,2);
        D = D.*(size(X,2)./sum(~isnan(X),2));
        D = sqrt(D);
    case 'cit'
        D = bsxfun(@minus,X,C(L,:));
        D = nansum(abs(D),2);
        D = D.*(size(X,2)./sum(~isnan(X),2));
    case 'sqcit'
        D = bsxfun(@minus,X,C(L,:));
        D = nansum(abs(D),2);
        D = D.*(size(X,2)./sum(~isnan(X),2));
        D = D.^2;
end
        
end

function D = nandistfuncp2_exp(C, Xi, L, distance, sx) 
% Description: 
% Computes distances between observations with missing values and nearest 
% centroids. Uses expected distance strategy for treating missing values. 
%
% Inputs:
%        C - Cluster centroids 
%       Xi - Imputed data set using conditional mean 
%        L - Cluster labels 
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
    sx = zeros(size(Xi));
end
switch distance
    case 'euc'
        % Compute omega
        omega = bsxfun(@minus,Xi,C(L,:));
        omega = sum(omega.^2,2) + sum(sx,2);  
        % Compute variance
        Ex = Xi;
        Ex2 = Xi.^2 + sx;
        Ex3 = Xi.^3 + 3*Xi.*sx;
        Ex4 = Xi.^4 + 6*(Xi.^2).*sx + 3*sx.^2;
        Y = C(L,:);
        Ey = Y;
        Ey2 = Y.^2;
        Ey3 = Y.^3;
        Ey4 = Y.^4;
        var = sum(Ex4 + Ey4 - 4*Ex3.*Ey - 4*Ex.*Ey3 + 6*Ex2.*Ey2,2) - ...
            sum((Ex2 - 2*Ex.*Ey + Ey2).^2,2);
        var(var<0.0000001) = 0; 
        % Compute EED
        m = (omega.^2)./var;
        D = exp(gammaln(m+0.5)-gammaln(m));
        D = D.*((omega./m).^(0.5));
        ind = isnan(D);
        D(ind) = sqrt(omega(ind));
    case 'sqe'
        % Compute ESD
        sx = sum(sx,2);
        D = bsxfun(@minus,Xi,C(L,:));
        D = sum(D.^2,2) + sx;
end

end

function D = nanmatrixdist_ads(X, Y, distance)
% Description: 
% Computes distances between two equal sized matrices using 
% available data strategy is used for treating missing values. 
%
% Inputs:
%         X - First data set
%         Y - Second data set
%  distance - Selected distance metric 
%             Alternatives: 
%             'sqe' - squared Euclidean distance
%             'euc' - Euclidean distance 
%             'cit' - City block distance 
%             'sqcit' - squared City block distance 
%
% Output:
%         D - Distance matrix
switch distance
    case 'sqe'
        D = bsxfun(@minus,X,Y);       
        D = nansum(D.^2,2);
    case 'euc'
        D = bsxfun(@minus,X,Y);       
        D = sqrt(nansum(D.^2,2));
    case 'cit'
        D = bsxfun(@minus,X,Y);      
        D = nansum(abs(D),2);
    case 'sqcit'
        D = bsxfun(@minus,X,Y);
        D = nansum(abs(D).^2,2);
end
             
end

function D = nanmatrixdist_pds(X, Y, distance)
% Description: 
% Computes distances between two equal sized matrices using 
% partial distance strategy for treating missing values.
%
% Input:
%         X - First data set 
%         Y - Second data set
%  distance - Selected distance metric 
%             Alternatives: 
%             'sqe' - squared Euclidean distance 
%             'euc' - Euclidean distance 
%             'cit' - City block distance  
%            'scit' - squared City block distance  
%
% Output:
%         D - Distance matrix
switch distance    
    case 'sqe'
        D = bsxfun(@minus,X,Y);
        m = sum(~isnan(D),2);        
        D = nansum(D.^2,2);
        D = D.*(size(X,2)./m);
    case 'euc'
        D = bsxfun(@minus,X,Y);
        m = sum(~isnan(D),2);        
        D = nansum(D.^2,2);
        D = sqrt(D.*(size(X,2)./m));
    case 'cit'
        D = bsxfun(@minus,X,Y);
        m = sum(~isnan(D),2);        
        D = nansum(abs(D),2);
        D = D.*(size(X,2)./m);
    case 'sqcit'
        D = bsxfun(@minus,X,Y);
        m = sum(~isnan(D),2);        
        D = nansum(abs(D),2);
        D = (D.*(size(X,2)./m)).^2;        
end
             
end

function D = nanmatrixdist_exp(X, Y, distance, sx, sy)
% Description: 
% Computes distances between two equal sized matrices using expected 
% distances for treating missing values. 
%
% Input:
%          X - First imputed data set
%          Y - Second imputed data set  
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
        % Compute omega
        s = sum(sx,2) + sum(sy,2);
        omega = bsxfun(@minus,X,Y);       
        omega = sum(omega.^2,2);
        omega = omega + s;
        % Compute variance
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
        % Compute EED
        m = (omega.^2)./var;
        D = exp(gammaln(m + 0.5) - gammaln(m));
        D = D.*((omega./m).^(0.5));
        ind = isnan(D);
        D(ind) = sqrt(omega(ind));
    case 'sqe'
        % Compute ESD
        sx = sum(sx,2); 
        sy = sum(sy,2);
        s = sx + sy;  
        D = bsxfun(@minus,X,Y);       
        D = sum(D.^2,2);
        D = D + s;            
end
             
end

function C = nanmean(X, weights, varargin)
% Description: 
% Computes mean value of input data set. 
%
% Input:
%         X - Input data set
%   weights - Weights of input data set  
%
% Output:
%         C - Mean value of data set  
%
nan = isnan(X);
X(nan) = 0;
C(1,:) = sum(weights.*X,1) ./ sum(weights.*~nan,1);

end

function C = nanmedian(X, weights, varargin)
% Description: 
% Computes median value of input data set. 
%
% Input:
%         X - Input data set 
%   weights - Weights of input data    
%
% Output:
%         C - Median value of data set
%
n = size(X,2);
wMat = weights.*~isnan(X);
nn = sum(wMat,1);
for j = 1:n
    [X(:,j), idx] = sort(X(:,j));
    wMat(:,j) = wMat(idx,j); 
end
wMat = wMat./sum(wMat,1);
wMat(isinf(wMat)) = 0;
C = nan(1,n);
for j = 1:n
    jj = find(cumsum(wMat(:,j),1)>=0.5,1);
    cumsumW = cumsum(wMat(1:jj,j));
    cumsumW = round(gather(cumsumW(end)),10);
    if isempty(jj)
        C(1,j) = NaN;
    elseif mod(nn(j),2)==0 && cumsumW==0.5
        C(1,j) = 0.5*(X(jj,j)+X(jj+1,j));
    else
        C(1,j) = X(jj,j);
    end
end

end

function C = nanspatialmedian_ads(X, weights, varargin)
% Description: 
% Computes spatial median value of input data set. 
% Uses available data strategy for treating missing values.
%
% Inputs: 
%         X - Input data set
%   weights - Weights of input data    
%
% Output: 
%         C - Spatial-median value of data set
%
u = zeros(1,size(X,2));
for i = 1:size(X,2)
    I = ~isnan(X(:,i));
    u(i) = median(X(I,i));
end
max_iter = 100;
tol = 1e-5;
P = ~isnan(X);
iters = 0;
w = 1.5;
X(isnan(X)) = 0;
while iters < max_iter
    iters = iters + 1;
    D = P.*(X-u);
    a = weights./sqrt(sum(D.^2,2)+sqrt(eps));
    a = P.*a;
    ax = sum(a.*X,1);
    v = (1./sum(a)).*ax;
    u1 = u + w*(v-u);
    if norm(u1-u,inf) < tol
        break;
    end
    u = u1;
end
C = u1;

end

function C = nanspatialmedian_pds(X, weights, varargin)
% Description: 
% Computes spatial median value of input data set. 
% Uses partial distance strategy for treating missing values.
%
% Inputs: 
%         X - Input data set
%   weights - Weights of input data   
%
% Output: 
%         C - Spatial-median value of data set
%
u = zeros(1,size(X,2));
for i = 1:size(X,2)
    I = ~isnan(X(:,i));
    u(i) = median(X(I,i));
end
max_iter = 100;
tol = 1e-5;
P = ~isnan(X);
iters = 0;
w = 1.5;
Xorg = X;
X(isnan(X)) = 0;
while iters < max_iter
    iters = iters + 1;
    D = P.*(X-u);
    alpha = size(X,2)./sum(~isnan(Xorg),2);
    a = weights./sqrt(alpha.*sum(D.^2,2)+sqrt(eps));
    a = P.*a;
    ax = sum(a.*X,1);
    v = (1./sum(a)).*ax;
    u1 = u + w*(v-u);
    if norm(u1-u,inf) < tol
        break;
    end
    u = u1;
end
C = u1;

end

function C = nanspatialmedian_exp(X, weights, Xi, sx)  
% Description: 
% Computes spatial median value of input data set. 
% Uses expected distance strategy for treating missing values.
%
% Inputs: 
%         X - Input data set with missing values
%   weights - Weights of input data   
%        Xi - Conditional mean imputed data set 
%        sx - Variances of imputed data vectors       
%
% Output: 
%         C - Spatial-median value of data set
%
u = zeros(1,size(X,2));
for i = 1:size(X,2)
    I = ~isnan(X(:,i));
    u(i) = median(X(I,i));
end
max_iter = 100;
tol = 1e-5;
P = ~isnan(X);
iters = 0;
w = 1.5;
X(isnan(X)) = 0;
while iters < max_iter
    iters = iters + 1;
    D = Xi-u;
    a = weights./sqrt(sum(D.^2,2)+sum(sx,2)+sqrt(eps));
    a = P.*a;
    ax = sum(a.*X,1);
    v = (1./sum(a)).*ax;
    u1 = u + w*(v-u);
    if norm(u1-u,inf) < tol
        break;
    end
    u = u1;
end
C = u1;

end


