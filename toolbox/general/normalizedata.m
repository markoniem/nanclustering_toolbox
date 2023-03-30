function X = normalizedata(X, varargin)
% Description: 
% Rescales columns of input data based on min-max or z-score 
% normalization.  
%
% Inputs:
%         X - Input data set
%  'Method' - Pre-processing method 
%             Alternatives: 
%             'min-max' - Min-max normalization 
%             'zscore' - Z-score normalization
%              Default is: min-max.  
%   'Range' - Range of min-max scaling. Default is: [-1, 1].  
%
% Output:
%         X - Pre-procesed data set
%
pnames = {'method' 'range'};
dflts =  {'min-max' [-1,1]};
[method, range] = internal.stats.parseArgs(pnames, dflts, varargin{:});
methodNames = {'min-max','zscore'};
if isnumeric(method) || ~ismember(method,methodNames)
    error('Invalid scaling method');
end
%
switch method
    case 'min-max'
        a = range(1);
        b = range(2);
        min_x = min(X);
        max_x = max(X);
        cofs = (b - a) ./ (max_x - min_x);
        if (b - a) <= sqrt(eps), error('Range must be non-negative!'); end
        X = bsxfun(@minus,X,min_x);
        X = a + bsxfun(@times,X,cofs);
    case 'zscore'
        X = bsxfun(@minus,X,mean(X,'omitnan'));
        X = bsxfun(@rdivide,X,std(X,'omitnan'));
    otherwise
        fprintf('Nothing to be done...\n')
end

end

