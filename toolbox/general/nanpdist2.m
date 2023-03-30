function D = nanpdist2(Xm, Ym, varargin)
% Description:
% Computes pairwise distances between two data sets. Function assumes that 
% data contains missing values. Available data, partial, expected squared 
% euclidean, and expected euclidean distances are available options. 
%
% Inputs:
%                  Xm - Input data set with missing values 
%                  Ym - Second data set with missing values 
%          'Distance' - Selected distance metric. Default is 'euc'.
%                       Alternatives: 
%                       'sqe' - squared Euclidean distance
%                       'euc' - Euclidean distance
%                       'cit' - City block distance 
%      'TreatMissing' - Used distance computation function. 
%                       Alternatives: 
%                       'ads' - Available data 
%                       'pds' - Partial distance 
%                       'esd' - Expected squared Euclidean
%                       'eed' - Expected Euclidean
%                       Default is: 'ads'
%
% Output:
%                   D - Computed pairwise distances
%
pnames = {'distance' 'treatmissing'};
dflts =  {'euc', 'ads'};
[distance, treatmissing] = internal.stats.parseArgs(pnames, dflts, varargin{:});
distNames = {'sqe','euc','cit'};
treatmissingNames = {'ads','pds','esd','eed'};
if isnumeric(distance) || ~ismember(distance,distNames)
    error('Invalid distance');
end
if isnumeric(distance) || ~ismember(treatmissing,treatmissingNames)
    error('Invalid treatmissing');
end
if ((strcmp(treatmissing,'esd')||strcmp(treatmissing,'eed')) && strcmp(distance,'cit'))
    error('Expected City block is not supported');
end
%
if strcmp(treatmissing,'ads')
    nandistfunc = @nandistfunc_ads;
    sx = [];
    sy = [];
elseif strcmp(treatmissing,'pds')
    nandistfunc = @nandistfunc_pds;
    sx = [];
    sy = [];
else
    nandistfunc = @nandistfunc_exp;
    [Xm, sx, success] = ecmnmlefunc(Xm);
    [Ym, sy, success2] = ecmnmlefunc(Ym);
    if ~success || ~success2 
        error('Error for computing parameters of conditional distribution');
    end
end
%
D = nandistfunc(Xm, Ym, distance, sx, sy);

end

function D = nandistfunc_ads(X, Y, distance, varargin)
% Description: 
% Computes distances between observations with missing values in two data 
% matrices. Uses available data strategy for treating missing values.
%
% Inputs:
%         X - First data set with missing values
%         Y - Second data set set with missing values
%  distance - Selected distance metric 
%             Alternatives: 
%             'sqe' - squared Euclidean distance
%             'cit' - City block distance 
%             'euc' - Euclidean distance 
%
% Output:
%         D - Distance matrix 
%
switch distance
    case 'sqe'
        I = isnan(Y);
        Y(I) = 0;
        I2 = isnan(X);
        X(I2) = 0;
        D = ((X.^2)*~I'-2*X*Y'+~I2*(Y.^2)');
        D(D<eps) = 0;
    case 'cit'
        [nx, p] = size(Y);
        [nc, ~] = size(X);
        D = zeros(nc,nx,'double');
        for i = 1:nc
            dsq = zeros(nx,1,'double');
            for q = 1:p
                dsq1 = abs(Y(:,q)-X(i,q));
                dsq1(isnan(dsq1)) = 0;
                dsq = dsq + dsq1;
            end
            D(i,:) = dsq;
        end
    case 'euc'
        I = isnan(Y);
        Y(I) = 0;
        I2 = isnan(X);
        X(I2) = 0;
        D = ((X.^2)*~I'-2*X*Y'+~I2*(Y.^2)');
        D(D<eps) = 0;
        D = sqrt(D);
end

end

function D = nandistfunc_pds(X, Y, distance, varargin)
% Description: 
% Computes distances between two matrices using partial
% distance strategy for treating missing values.
%
% Inputs: 
%        X - First data set with missing values
%        Y - Second data set with missing values
%  distance - Selected distance metric 
%             Alternatives: 
%             'sqe' - squared Euclidean distance      
%             'cit' - City block distance  
%             'euc' - Euclidean distance 
%             
% Output:
%        D - Distance matrix 
switch distance
    case 'sqe'
        p = size(X,2);
        I = isnan(Y);
        Y(I) = 0;
        I2 = isnan(X);
        X(I2) = 0;
        D = ((X.^2)*~I'-2*X*Y'+~I2*(Y.^2)');
        D(D<eps) = 0;
        D = D.*(p./(double(~I2)*double(~I')));
    case 'cit'
        [nx, p] = size(Y);
        [nc, ~] = size(X);
        D = zeros(nc,nx,'double');
        for i = 1:nc
            dsq = zeros(nx,1,'double');
            m = p*ones(nx,1);
            for q = 1:p
                dsq1 = abs(Y(:,q)-X(i,q));
                isn = isnan(dsq1);
                m(isn) = m(isn) - 1;
                dsq1(isn) = 0;
                dsq = dsq + dsq1;
            end
            D(i,:) = dsq.*(p./m);
        end
    case 'euc'
        p = size(X,2);
        I = isnan(Y);
        Y(I) = 0;
        I2 = isnan(X);
        X(I2) = 0;
        D = ((X.^2)*~I'-2*X*Y'+~I2*(Y.^2)');
        D(D<eps) = 0;
        D = D.*(p./(double(~I2)*double(~I')));
        D = sqrt(D);
end

end

function D = nandistfunc_exp(X, Y, distance, sx, sy)
% Description: 
% Compute distances between two matrices. Uses expected distance 
% estimation for treating missing values. 
%
% Inputs: 
%         X - Imputed data set using conditional mean
%         Y - Second imputed data set using conditional mean
%  distance - Selected distance metric 
%             Alternatives: 
%             'euc' - Euclidean distance 
%             'sqe' - squared Euclidean distance
%
% Output:
%       D - Distance matrix
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
        omega = pdist2(X,Y).^2;
        omega = bsxfun(@plus,omega,sum(sx,2));
        omega = bsxfun(@plus,omega,sum(sy,2)');  
        % Compute variance
        Ex = X;
        Ex2 = X.^2 + sx;
        Ex3 = X.^3 + 3*X.*sx;
        Ex4 = X.^4 + 6*(X.^2).*sx + 3*sx.^2;
        Ey = Y;
        Ey2 = Y.^2 + sy;
        Ey3 = Y.^3 + 3*Y.*sy;
        Ey4 = Y.^4 + 6*(Y.^2).*sy + 3*sy.^2;
        I = ones(size(Y));
        I2 = ones(size(X));
        var = Ex4*I' + I2*Ey4' - 4*Ex3*Ey' - 4*Ex*Ey3' + 6*Ex2*Ey2' - ...
            (Ex2*I' - 2*Ex*Ey' + I2*Ey2').^2;
        var(var<0.0000001) = 0;
        % Compute EED
        m = (omega.^2)./var;
        D = exp(gammaln(m+0.5) - gammaln(m));
        D = D.*((omega./m).^(0.5));
        ind = isnan(D);
        D(ind) = sqrt(omega(ind));  
    case 'sqe'
        % Compute ESD
        D = pdist2(X,Y).^2;
        sx = sum(sx,2);
        sy = sum(sy,2);
        D = bsxfun(@plus,D,sx);
        D = bsxfun(@plus,D,sy');
end
        
end

