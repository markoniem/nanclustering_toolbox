function plot_indices(clusterInd,indices_values)
% Description: 
% Plots curves of cluster validation indices   
%
% Inputs:
%       clusterInd - Selected cluster validation indices
%       indices_values - Values of cluster validation indices
%
for i = 1:size(indices_values,1)
    figure;
    plot(indices_values(i,:));
    title(char(clusterInd{i,1}));
end

end

