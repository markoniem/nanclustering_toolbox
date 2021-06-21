function D = PDS(X)
% Description: 
% Compute pairwise distance between observations with missing values by 
% using partial distance strategy.
%
% Function call:
%         D = PDS(X)
%
% Input:
%         X - Input data set 
%
% Output:
%         D - Calculated distances
%
n = size(X,1);
p = size(X,2);
Y = X;
I = isnan(Y);
Y(I) = 0;
I2 = isnan(X);
X(I2) = 0;
D = ((X.^2)*~I'-2*X*Y'+~I2*(Y.^2)');
D(D<eps) = 0;
D = D.*(p./(double(~I2)*double(~I')));
B = tril(ones(n,n),-1);
D = sqrt(D(B==1));

end

