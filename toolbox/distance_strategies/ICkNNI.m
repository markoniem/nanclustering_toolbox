function D = ICkNNI(X)
% Description: 
% Compute pairwise distance between observations with missing values using  
% Incomplete-Case kNN Imputation. 
%
% Function call:
%         D = ICkNNI(X)
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
[N, n] = size(X);
Im = find(any(isnan(X),2));
Ximp = X;
for i = 1:N
    if ~ismember(i,Im)
        continue;
    end
    Xmi = X(i,:);
    miss = find(isnan(Xmi));
    notmiss = find(~isnan(Xmi));
    D = nandistfunc(Xmi,X,distance);
    [~, idx] = sort(D);
    neighbors = zeros(length(miss),K);
    Kcounter = zeros(length(miss),1);
    for j = 2:length(idx)
        Xj = X(idx(j),:);
        miss2 = find(isnan(Xj));
        if sum(ismember(miss2,notmiss)) > 0
            continue;
        else
            for k = 1:length(miss)
                if ismember(miss(k),miss2) || Kcounter(k) == K
                    continue;
                else
                    Kcounter(k) = Kcounter(k) + 1;
                    neighbors(k, Kcounter(k)) = idx(j);
                end
            end
        end
        if isempty(find(neighbors==0,1))
            break;
        end
    end
    for j = 1:length(miss)
        if Kcounter(j) > 0
            neighbors1 = neighbors(j,1:Kcounter(j));
            Xmi(miss(j)) = sum(X(neighbors1,miss(j)))/Kcounter(j);
        end
    end
    Ximp(i,:) = Xmi;
end
if ~isempty(isnan(Ximp))
    for i = 1:n
        X1 = Ximp(:,i);
        X1(isnan(X1)) = mean(X(:,i),'omitnan');
        Ximp(:,i) = X1;
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

