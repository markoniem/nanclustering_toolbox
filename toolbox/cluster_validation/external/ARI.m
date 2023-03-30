function ARIval = ARI(x,y)
% Description: 
% Computes adjusted Rand index value.
%
% Inputs:
%         x - Real labels of clusters
%         y - Predicted labels of clusters  
% Output:
%    ARIval - Value of adjusted Rand index  
%
N = length(x);
n = zeros(max(x),max(y));
for i = 1:max(x)
    for j = 1:max(y)
        G1 = find(x==i);
        G2 = find(y==j);
        n(i,j) = length(intersect(G1,G2));
    end
end
index = 0;
a = 0;
b = 0;
for i=1:max(x)
    for j=1:max(y)
        if n(i,j) > 1
            index = index + nchoosek(n(i,j),2);
        end
    end
end
for i = 1:max(x)
    if sum(n(:,i)) > 1
        a = a + nchoosek(sum(n(:,i)),2);
    end
end
for i = 1:max(y)
    if sum(n(i,:)) > 1
        b = b + nchoosek(sum(n(i,:)),2);
    end
end
NN = index - a*b/nchoosek(N,2);
DD = (a + b)/2 - a*b/nchoosek(N,2);
if (NN == 0 && DD==0)
    ARIval = 1;
else
    ARIval = NN/DD;
end

end

