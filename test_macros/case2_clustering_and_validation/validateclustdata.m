% Description: 
% Perform cluster validation and saves cluster validation results to file. 
%
% Reference results are available in the following article: 
% M. Niemelä and T. Kärkkäinen, "Improving clustering and cluster validation 
%     with missing data using distance estimation methods," In: Computational 
%     Sciences and Artificial Intelligence in Industry: New digital technologies 
%     for solving future societal and economical challenges, Intelligent Systems, 
%     Control and Automation: Science and Engineering. Springer Verlag, 2021. 
%     (12 pages, to appear).
%
% Inputs (see selected values from 'parameters.m' file):
%        datasets - Input data set(s)
%    datamatrices - Clustered data sets which consisted of predefined numbers 
%                   of missing values (loaded from file)
%         centers - Obtained cluster centroids (loaded from file)
%          labels - Obtained cluster labels for each observation (loaded from file)
%        probmiss - Probabilities of missing values generated to data set(s)
%        distance - Selected distance metric. 
%                   Alternatives: 
%                   'sqe' - squared Euclidean distance
%                   'euc' - Euclidean distance
%                   'cit' - City block distance 
%     distFuncInd - Used distance estimation strategy for indices 
%                   Alternatives:
%                   'ads' - Available data 
%                   'pds' - Partial distance 
%                   'exp' - Expected distance (ESD/EED)
%      clusterInd - Selected cluster validation indices
%                   Alternatives (multiple indices can be selected): 
%                   @CalinskiHarabasz - Calinski-Harabasz 
%                   @DaviesBouldin - Davies-Bouldin 
%                   @DaviesBouldin2 - Davies-Bouldin*
%                   @GenDunn - Generalized Dunn
%                   @kCE - kCE-index
%                   @PBM - Pakhira-Bandyopadhyay-Maulik
%                   @RayTuri - Ray-Turi
%                   @Silhouette - Silhouette
%                   @Silhouette2 - Silhouette*
%                   @WB - WB-index
%                   @WemmertGancarski - Wemmert-Gancarski
%                   Note: Used distance metric(s) for the selected indices can 
%                   be optionally separately specified. 
%                   Use then the following format typing: 
%                   [{@index1, distance1}; {@index2, distance2}; ...], e.g.,
%                   [{@CalinskiHarabasz, 'sqe'}; {@DaviesBouldin, 'euc'}]
%                   In general case, if specific distances were not given, 
%                   the 'distance' parameter will be used (see definition above). 
%
% Outputs:
%
%  indices_values - Calculated index values for each number of clusters. 
%                   Indices will be saved to four dimensional matrix 
%                   such as first dimension is for cluster validation indices, 
%                   second dimension is for different number of clusters 
%                   (starting from 2), third dimension correspond different 
%                   data sets, and fourth dimension is for different number 
%                   of missing values.
% DyymmddThhm.xlsx - Recommended number of clusters by cluster validation 
%                    indices saved to .xlsx file format. File is named based 
%                    on execution time of macro. 
load('../../benchmark_data/CSAI_data/datamatrices');
load('centers');
load('labels');
addpath('../../toolbox/cluster_indices');
params = params();
datasets = params.datasets;
probmiss = params.probmiss;
distance = params.distance;
distFuncInd = params.distFuncInd;
clusterInd = params.clusterInd;
indices_values = zeros(length(clusterInd),size(centers,1)+1, ... 
                    length(datasets),length(probmiss)); 
for i = 1:length(datasets)
    fprintf('Data set: %s\n',datasets{i})
    for j = 1:length(probmiss)
        fprintf('Missing values: %.2f %% \n',probmiss(j));  
        Xm = datamatrices{i,j};
        indices_values(:,:,i,j) = cluster_validation(Xm,centers(:,i,j), ... 
                 labels(:,i,j),distance,clusterInd,distFuncInd);
    end
    fprintf('=================\n');
end
save('indices_values','indices_values');
%
% Save results to .xlsx file
%
[~, indices_results]= min(indices_values,[],2);
results = zeros(length(clusterInd),length(datasets),length(probmiss));
for i = 1:length(datasets), for j = 1:length(probmiss), results(:,i,j) = ...
            indices_results(:,:,i,j); end, end
for i = 1:length(datasets), datasets{i} = sprintf('Data set: %s',datasets{i}); end
clusterIndchar = cell(size(clusterInd,1),1);
for i = 1:length(clusterInd), clusterIndchar{i} = sprintf('%s',char(clusterInd{i})); end
for i = 1:size(results,3)    
    R = results(:,:,i);
    T = table(clusterIndchar,R);
    resultsTable = splitvars(T);
    resultsTable.Properties.VariableNames = [' ', datasets'];
    writetable(resultsTable,sprintf('D%sT%s.xlsx', ...
                datestr(now,'yymmdd'), datestr(now,'HHMM')), 'Sheet', ... 
                        sprintf('Missing values %2d%%',100*probmiss(i)));
end

