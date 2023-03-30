function ACCval = ACC(x,y)
% Description: 
% Computes value of clustering accuracy
%
% Inputs:
%         x - Real labels of clusters
%         y - Predicted labels of clusters  
% Output:
%    ACCval - Value of accuracy
%
I = find(x==y);
ACCval = length(I)/length(x); 

end

