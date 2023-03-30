% Description: 
% Performs clustering and cluster validation using high-dimensional data. 
% Results are available in a xlsx file. 
%
% Inputs: 
%       See parameters function
%
% Outputs:
%         centers - Cluster centers. Centers will be saved to three 
%                   dimensional cell array. First dimension is for different 
%                   cluster numbers, second dimension is for different data 
%                   sets, and third dimension is for different percentages 
%                   of missing values in data.
%          labels - Cluster labels for each observation. Labels will be saved 
%                   to three dimensional cell array. First dimension is for 
%                   different cluster numbers, second dimension is for 
%                   different data sets, and third dimension is 
%                   for different percentages of missing values in data.
%  indices_values - Computed index values for each number of clusters. 
%                   Indices will be saved to four dimensional matrix. 
%                   First dimension is for cluster validation indices, 
%                   second dimension is for different numbers of clusters 
%                   (starting from 2), third dimension is for different 
%                   data sets, and fourth dimension is for different number 
%                   of missing values.
%
clear, clc, close all;
addpath('../datasets/multidim');
addpath('../toolbox/general');
addpath('../toolbox/clustering');
addpath('../toolbox/cluster_validation');
params = parameters();
datasets = params.datasets;  
correctResults = params.correctResults;               
probmiss = params.probmiss;
maxClusters = params.maxClusters;
treatMissing = params.treatMissing;
treatMissingValidation = params.treatMissingValidation;
distance = params.distance;
replicates = params.replicates;
initcrit = params.initcrit;
useprevCent = params.useprevCent;
pipelined = params.pipelined;
showProgression = params.showProgression;
clusterInd = params.clusterInd;
%
centers = cell(maxClusters-1,length(datasets),length(probmiss));
labels = cell(maxClusters-1,length(datasets),length(probmiss));
indices_values = zeros(length(clusterInd),maxClusters, ... 
                    length(datasets),length(probmiss)); 
for i = 1:length(datasets)
    load(datasets{i});
    fprintf('Data set: %s\n',datasets{i})
    for j = 1:length(probmiss)
        fprintf('Missing values: %.2f %% \n',probmiss(j));  
        X = normalizedata(X,'Method','min-max','Range',[-1, 1]);
        [Xm, I] = genmissdata(X,probmiss(j)); 
        Xm = Xm(I,:);
        % Some complete observations are required for the initialization 
        % step of the clustering 
        p = uint64(0.01*size(Xm,1));
        I = ~any(isnan(Xm),2);
        if sum(I) < p
            I = randi(size(Xm,1),p,1);
            Xm(I,:) = X(I,:);
        end
        [centers(:,i,j), labels(:,i,j)] = iterative_kcentroids(Xm,'MaxClusters',maxClusters, ...
            'TreatMissing',treatMissing,'Distance',distance,'Replicates',replicates,'Initcrit', ...
            initcrit, 'Pipelined', pipelined, 'UsePrevCent', useprevCent, 'ShowProgression', showProgression);
        indices_values(:,:,i,j) = cluster_validation(Xm,centers(:,i,j),labels(:,i,j), ...
            'TreatMissing',treatMissingValidation,'Distance',distance,'ClusterInd',clusterInd);
    end
    fprintf('=================\n');
end
save('centers','centers');
save('labels','labels');
save('indices_values','indices_values');
%
% Save results to .xlsx file
%
[~, indices_results]= min(indices_values,[],2);
correctNumbers = zeros(length(clusterInd),length(probmiss));
for k = 1:length(probmiss)
    for i = 1:length(clusterInd)
        for j = 1:length(correctResults)
            if indices_results(i,1,j,k) == correctResults(j)
                correctNumbers(i,k) = correctNumbers(i,k) + 1;
            end
        end
    end
end
results = zeros(size(clusterInd,1),length(datasets),length(probmiss));
for i = 1:length(datasets), for j = 1:length(probmiss), results(:,i,j) = ...
            indices_results(:,:,i,j); end, end
for i = 1:length(datasets), datasets{i} = sprintf('Data set: %s',datasets{i}); end
datasets{length(datasets)+1} = 'Correct';
clusterIndchar = cell(size(clusterInd,1),1);
for i = 1:size(clusterInd,1), clusterIndchar{i} = sprintf('%s',char(clusterInd{i,1})); end
for i = 1:size(results,3)    
    R = results(:,:,i);
    R = [R, correctNumbers(:,i)];
    T = table(clusterIndchar,R);
    resultsTable = splitvars(T);
    resultsTable.Properties.VariableNames = [' ', datasets'];
    writetable(resultsTable,'partC_Results.xlsx', 'Sheet', ... 
                        sprintf('Missing values %2d%%',100*probmiss(i)));
end
restoredefaultpath;

function params = parameters()
% Description: 
% Return initial parameters for test macro. 
%
% Output:
%     params - Parameters of test macro
%
params.datasets = {'K15_Nk400_M10_dc0.9.mat'; 'K15_Nk400_M10_dc0.8.mat';
                   'K15_Nk400_M10_dc0.7.mat'; 'K15_Nk400_M10_dc0.6.mat';
                   'K15_Nk400_M50_dc0.9.mat'; 'K15_Nk400_M50_dc0.8.mat';
                   'K15_Nk400_M50_dc0.7.mat'; 'K15_Nk400_M50_dc0.6.mat';
                   'K15_Nk400_M100_dc0.9.mat'; 'K15_Nk400_M100_dc0.8.mat';
                   'K15_Nk400_M100_dc0.7.mat'; 'K15_Nk400_M100_dc0.6.mat'};                 
params.correctResults = [15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15];
params.probmiss = [0.00; 0.05; 0.10; 0.20];
params.maxClusters = 20;
params.treatMissing = 'ads';
params.treatMissingValidation = 'ads';
params.distance = 'euc';
params.replicates = 100;
params.initcrit = 'kmeans++';
params.useprevCent = true;
params.pipelined = true;
params.showProgression = true;
params.clusterInd = {@CalinskiHarabasz; @DaviesBouldin; @DaviesBouldin2; ...
    @GenDunn; @kCE; @PBM; @RayTuri; @Silhouette; @WB; @WemmertGancarski};
% params.clusterInd = [{@CalinskiHarabasz, 'sqe'};{@DaviesBouldin, 'euc'}; ...
%     {@DaviesBouldin2, 'euc'}; {@GenDunn, 'euc'}; {@kCE, 'sqe'}; ...
%     {@PBM, 'euc'}; {@RayTuri, 'sqe'}; {@Silhouette, 'euc'}; ...
%     {@WB, 'sqe'}; {@WemmertGancarski, 'euc'}];

end
