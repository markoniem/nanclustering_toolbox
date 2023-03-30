function kCEVal = kCE(X, C, L, weights, Xi, sx, varargin)
% Description: 
% Computes kCE-index value. Three options are available 
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
%        kCEVal - Value of kCE-index 
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
        nandistfuncp2 = @nandistfuncp2_ads;
    elseif strcmp(treatmissing,'pds')
        nandistfuncp2 = @nandistfuncp2_pds;
    end    
    k = size(C,1);
    switch distance
        case 'euc'
            kCEVal = k*nansum(weights.*nandistfuncp2(C,X,L,'sqe'),1);
        case 'sqe'
            kCEVal = k*nansum(weights.*nandistfuncp2(C,X,L,distance),1);
        case 'cit'
            kCEVal = k*nansum(weights.*nandistfuncp2(C,X,L,'sqcit'),1);
    end
else
    k = size(C,1);
    kCEVal = k*nansum(weights.*nandistfuncp2_exp(C,Xi,L,'sqe',sx),1);
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


