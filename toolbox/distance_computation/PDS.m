function [D, success] = PDS(X,varargin)
% Description: 
% Compute pairwise distance between observations with missing values using  
% partial distance strategy. 
%
% Input:
%          X - Input data set 
% 'Distance' - Selected distance metric. Default is 'euc'.
%              Alternatives: 
%              'sqe' - squared Euclidean distance
%              'euc' - Euclidean distance
%              'cit' - City block distance 
%
% Output:
%         D - Calculated distances
%   success - Indicator flag for successfully finished algorithm  
%
success = 1;
pnames = {'distance'};
dflts =  {'euc'};
[distance] = internal.stats.parseArgs(pnames, dflts, varargin{:});
distNames = {'sqe','euc','cit'};
if isnumeric(distance) || ~ismember(distance,distNames)
    error('Invalid distance');
end
D = nandistfunc(X, X, distance);
N = size(X,1);
B = tril(ones(N,N),-1);
D = D(B==1);

end

function D = nandistfunc(X, Y, distance)
% Description: 
% Computes distances between observations with missing values in two data 
% matrices. Uses partial strategy for treating missing values.
%
% Inputs:
%         X - First data set
%         Y - Second data set
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

