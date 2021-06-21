function D = kNNI(X)
% Description: 
% Compute pairwise distance between observations with missing values using  
% traditional kNN Imputation. 
%
% Function call:
%         D = kNNI(X)
%
% Input:
%         X - Input data set 
%
% Output:
%         D - Calculated distances
%
K = 5;
distance = 'euc';
Xorg = X;
I = ~all(isnan(X),2);
X = X(I,:);
Ximp = X;
Im = find(any(isnan(X),2));
for i = 1:length(Im)
    Xmi = X(Im(i),:);
    D = nandistfunc(Xmi,X,distance);
    [~, idx] = sort(D);
    cnt = 0;
    % In the case data vector consists missing values after
    % imputation, K will be increased one by one until vector is
    % fully complete.
    while sum(isnan(Xmi)) > 0
        knearests = X(idx(2:K+1+cnt),:);
        centroid = mean(knearests,'omitnan');
        ind = find(isnan(Xmi));
        Xmi(ind) = centroid(ind);
        Ximp(Im(i),:) = Xmi;
        cnt = cnt + 1;
    end
end
%
% Impute fully incomplete vectors
if (sum(I)<size(Xorg,1))
    Xcent = mean(Xorg,'omitnan');
    X0 = zeros(size(Xorg));
    cnt = 1;
    for i = 1:length(I)
        if I(i) == 1
            X0(i,:) = Ximp(cnt,:);
            cnt = cnt + 1;
        else
            X0(i,:) = Xcent;
        end
    end
    Ximp = X0;
end
D = pdist(Ximp)';

end

