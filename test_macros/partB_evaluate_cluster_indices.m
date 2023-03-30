% Description: 
% Evaluate performance of clustering and cluster validation algorithms with missing data. 
% Results are available in a xlsx file.  
%
% Reference results are available in the following article: 	
%  M. Niemelä and T. Kärkkäinen, "Improving clustering and cluster validation with 
%    missing data using distance estimation methods," In T. Tuovinen, J. Periaux, and 
%    P. Neittaanmäki, editors, Computational Sciences and Artificial Intelligence in 
%    Industry: New Digital Technologies for Solving Future Societal and Economical 
%    Challenges, pages 123–133. Springer International Publishing, 2022.
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
addpath('../datasets');
addpath('../toolbox/general');
addpath('../toolbox/clustering');
addpath('../toolbox/cluster_validation');
load('csai');
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
    fprintf('Data set: %s\n',datasets{i})
    for j = 1:length(probmiss)
        fprintf('Missing values: %.2f %% \n',probmiss(j));
        Xm = datamatrices{i,j};
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
    writetable(resultsTable,'partB_Results2.xlsx', 'Sheet', ... 
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
params.datasets = {'S1'; 'S2'; 'S3'; 'S4'; 'S5D2'; 'S2D2'; 'O200'; 'O2000'};
params.correctResults = [15, 15, 15, 15, 5, 2, 5, 5];
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

end
