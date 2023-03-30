function Silh = Silhouette(X, C, L, weights, Xi, sx, varargin)
% Description: 
% Computes Silhouette index value. Three options are available 
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
%        Silh - Value of Silhouette index 
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
if ~strcmp(treatmissing,'exp')
    if strcmp(treatmissing,'ads')
        nanmatrixdist = @nanmatrixdist_ads;
    elseif strcmp(treatmissing,'pds')
        nanmatrixdist = @nanmatrixdist_pds;    
    end    
    cnames = unique(L);
    k = length(cnames);
    if k == 1, Silh = Inf; return; end
    n = length(L);
    mbrs = (repmat(1:k,n,1) == repmat(L,1,k));
    avgDWithin = Inf(n,1);
    avgDBetween = Inf(n,k); 
    for j = 1:n
        distj = nansum(nanmatrixdist(X,X(j,:),distance),2);
        distj(distj == 0) = NaN;
        for i = 1:k
            if i == L(j)
                mbrs1 = mbrs;
                mbrs1(j,i) = 0;
                avgDWithin(j) = nanmean(distj(mbrs1(:,i)),weights(mbrs1(:,i)));
            else
                avgDBetween(j,i) = nanmean(distj(mbrs(:,i)),weights(mbrs(:,i)));
            end
        end
    end    
    minavgDBetween = min(avgDBetween, [], 2);
    Silh = (minavgDBetween - avgDWithin) ./ max(avgDWithin,minavgDBetween);
    Silh = nanmean(Silh,ones(size(Silh,1),1));
    % Inverse
    Silh = 1 / Silh;
else
    cnames = unique(L);
    k = length(cnames);
    if k == 1, Silh = Inf; return; end
    n = length(L);
    mbrs = (repmat(1:k,n,1) == repmat(L,1,k));
    avgDWithin = Inf(n,1);
    avgDBetween = Inf(n,k);
    Ex = Xi;
    Ex2 = Xi.^2 + sx;
    Ex3 = Xi.^3 + 3*Xi.*sx;
    Ex4 = Xi.^4 + 6*(Xi.^2).*sx + 3*sx.^2;
    for j = 1:n
        distj = nansum(nanmatrixdist2_exp(Xi,sx,Ex,Ex2,Ex3,Ex4,j),2);
        distj(distj == 0) = NaN;
        for i = 1:k
            if i == L(j)
                mbrs1 = mbrs;
                mbrs1(j,i) = 0;
                avgDWithin(j) = nanmean(distj(mbrs1(:,i)),weights(mbrs1(:,i)));
            else
                avgDBetween(j,i) = nanmean(distj(mbrs(:,i)),weights(mbrs(:,i)));
            end
        end
    end
    minavgDBetween = min(avgDBetween, [], 2);
    Silh = (minavgDBetween - avgDWithin) ./ max(avgDWithin,minavgDBetween);
    Silh = nanmean(Silh,ones(size(Silh,1),1));
    % Inverse
    Silh = 1 / Silh;    
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

function D = nanmatrixdist2_exp(X, sx, Ex, Ex2, Ex3, Ex4, idx)
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

