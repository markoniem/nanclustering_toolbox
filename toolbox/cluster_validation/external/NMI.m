function NMIval = NMI(x,y)
% Description: 
% Computes normalized mutual information value.
%
% Inputs:
%         x - Real labels of clusters
%         y - Predicted labels of clusters  
% Output:
%    NMIval - Value of normalized mutual information  
%
N = length(x);
Mx = zeros(N,max(x));
My = zeros(N,max(y));
for i =1:N
    Mx(i,x(i)) = 1;
    My(i,y(i)) = 1;
end
Pxy = nonzeros(Mx'*My/N);
Hxy = -dot(Pxy,log2(Pxy));
Px = nonzeros(mean(Mx,1));
Py = nonzeros(mean(My,1));
Hx = -dot(Px,log2(Px));
Hy = -dot(Py,log2(Py));
% Mutual information
MI = Hx + Hy - Hxy;
% Normalized mutual information
NMIval = MI/((Hx+Hy)/2);
NMIval = max(0,NMIval);

end

