% Description: 
% Performs K-centroids clustering using key points selection for selecting
% initial points of clustering.  The point are selected based on relatively 
% high density value (based on four nearest neighbors of data points) 
% and high density-based distance.
%
% Inputs: 
%       See parameters function
%
% Outputs:
% 
% new_datamatrices - Data sets after missing values generation and ICkNN 
%                    imputation. First dimension of 2D array correspond 
%                    different data sets and second dimension correspond 
%                    different degrees of missing values used in missing
%                    values generation.
%          centers - Cluster centers. Centers will be saved to three 
%                    dimensional cell array. First dimension is for different 
%                    cluster numbers, second dimension is for different data 
%                    sets, and third dimension is for different percentages 
%                    of missing values in data.
%        keypoints - Key points that are used in initialization of
%                    clustering. Data structure is same as the structure of
%                    the centers.
%           labels - Cluster labels for each observation. Labels will be saved 
%                    to three dimensional cell array. First dimension is for 
%                    different cluster numbers, second dimension is for 
%                    different data sets, and third dimension is 
%                    for different percentages of missing values in data.
%   indices_values - Computed index values for each number of clusters. 
%                    Indices will be saved to four dimensional matrix. 
%                    First dimension is for cluster validation indices, 
%                    second dimension is for different numbers of clusters 
%                    (starting from 2), third dimension is for different 
%                    data sets, and fourth dimension is for different number 
%                    of missing values.
clear, clc, close all;
addpath('../datasets');
addpath('../toolbox/general');
addpath('../toolbox/clustering');
addpath('../toolbox/cluster_validation');
load('csai');
%
params = parameters();
datasets = params.datasets;
correctResults = params.correctResults;
correctResults2 = params.correctResults2;
probmiss = params.probmiss;
maxClusters = params.maxClusters;
treatMissingValidation = params.treatMissingValidation;
distance = params.distance;
neighborsNum = params.neighborsNum;
showProgression = params.showProgression;
clusterInd = params.clusterInd;
%
new_datamatrices = cell(length(datasets),length(probmiss));
centers = cell(maxClusters-1,length(datasets),length(probmiss));
labels = cell(maxClusters-1,length(datasets),length(probmiss));
keypoints = cell(maxClusters-1,length(datasets),length(probmiss));
indices_values = zeros(length(clusterInd),maxClusters,length(datasets),length(probmiss)); 
%             
for i = 1:length(datasets)
    fprintf('Data set: %s\n',datasets{i})
    for j = 1:length(probmiss)
        fprintf('Missing values: %.2f %% \n',probmiss(j));  
        if strcmp(datasets{i},'iris_2D')||strcmp(datasets{i},'ecoli_2D')||strcmp(datasets{i},'seeds_2D')
            load(datasets{i});
            X = normalizedata(X,'Method','min-max','Range',[-1,1]);
            [Xm, I] = genmissdata(X,probmiss(j));
            Xm = Xm(I,:);
        else
            Xm = datamatrices{i,j};
        end
        if sum(isnan(Xm(:))) > 0
            Ximp = ICknnimpute(Xm, 'K', neighborsNum, 'Distance', distance);
            Xnew = Ximp;
        else
            Xnew = Xm;
        end
        new_datamatrices{i,j} = Xnew;
        [centers(:,i,j), labels(:,i,j), keypoints(:,i,j)] = iterative_keypoint_selection(Xnew, ...
                'MaxClusters',maxClusters,'Distance',distance,'ShowProgression',showProgression);
        indices_values(:,:,i,j) = cluster_validation(Xnew,centers(:,i,j),labels(:,i,j), ...
            'TreatMissing',treatMissingValidation,'Distance',distance,'ClusterInd',clusterInd);
    end
    fprintf('=================\n');
end
save('new_datamatrices','new_datamatrices');
save('centers','centers');
save('labels','labels');
save('keypoints','keypoints');
save('indices_values','indices_values');
%
% Save results to .xlsx file
%
[~, indices_results]= min(indices_values,[],2);
correctNumbers = zeros(length(clusterInd),length(probmiss));
for k = 1:length(probmiss)
    for i = 1:length(clusterInd)
        for j = 1:length(correctResults)
            if (indices_results(i,1,j,k) == correctResults(j) || indices_results(i,1,j,k) == correctResults2(j))
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
    writetable(resultsTable,'partB_Results4.xlsx', 'Sheet', ... 
                        sprintf('Missing values %2d%%',100*probmiss(i)));
end
restoredefaultpath;


function [centers, labels, keypoints] = iterative_keypoint_selection(X, varargin)
% Description: 
% Performs iteratively K-centroids clustering. Key points selection is used 
% for selecting initial points of clustering.
%
% Inputs:
%                  X - Input data set
%      'MaxClusters' - Maximum number of clusters. Default value is 20. 
%         'Distance' - Selected distance metric. Default value is 'euc'.
%                      Alternatives: 
%                      'sqe' - squared Euclidean distance
%                      'euc' - Euclidean distance  
%  'ShowProgression' - Indicator flag for presenting progression of clustering. 
%                      Default value is true.
% Outputs:
%            centers - Obtained cluster centers
%             labels - Obtained cluster labels
%
pnames = {'maxclusters' 'distance' 'showprogression'};
dflts =  {20 'euc' true};
[maxclusters, distance, showprogression] = internal.stats.parseArgs(pnames, dflts, varargin{:});
if (~isnumeric(maxclusters) || maxclusters < 2) 
    error('Invalid maxclusters');
end
%
centers = cell(maxclusters-1,1);
labels = cell(maxclusters-1,1);
keypoints = cell(maxclusters-1,1);
fprintf('Performing clustering...\n');
for k = 2:maxclusters
    if showprogression
        if k~= maxclusters, fprintf('k: %d, ',k); else, fprintf('k: %d\n',k); end 
        if mod(k,10) == 0 && k~= maxclusters, fprintf('\n'); end
    end
    kp1 = keypointsComp(X, distance, k, uint32(1.25*maxclusters));
    kp2 = keypointsComp2(X, distance, k, uint32(1.25*maxclusters));
    [~, ~, sumd1] = kcentroids(X, k, 'Distance', distance, 'Replicates', 100, 'Start', X(kp1,:));
    [~, ~, sumd2] = kcentroids(X, k, 'Distance', distance, 'Replicates', 100, 'Start', X(kp2,:));
    [~, I] = sort([sumd1; sumd2]);
    kp = [kp1, kp2];        
    kp = kp(:,I(1));
    keypoints{k-1} = X(kp,:);
    [L, C, ~] = kcentroids(X, k, 'Distance', distance, 'Replicates', 100, 'Start', X(kp,:));
    centers{k-1} = C;
    labels{k-1} = L;
end
fprintf('Done! \n');

end

function kp_idx = keypointsComp(X, dist, clusterNum, keyPointsNum)
% Description:
% Algorithm selects initial key points and removes points one-by-one 
% until clusterNum is reached. The point are removed based on 
% relatively low density and low density-based distance. Algorithm uses 
% median value for updating location of one key point during each iteration. 
%
% Inputs:
%            X - Input data set. Note: data must be complete.
%         dist - Selected distance metric. Default is: 'euc'.
%                Alternatives:
%                'sqe' - squared Euclidean distance
%                'euc' - Euclidean distance
%   clusterNum - The final number of clusters
% keyPointsNum - Initial number of key points
%
% Output:
%       kp_idx - Key point indexes
%
%Density
density = zeros(size(X,1),1);
for i = 1:size(X,1)
    Xi = X(i,:);
    D = pdist2(Xi,X,dist);
    [~, idx] = sort(D);
    knearests = X(idx(2:5),:);
    for j = 1:size(knearests,1)
        density(i) = density(i) + pdist2(Xi,knearests(j,:),dist);
    end
    density(i) = 1/density(i);
end
%
%Key points
sigma = inf(size(X,1),1);
kp = inf(size(X,1),1);
for i = 1:size(X,1)
    Xi = X(i,:);
    D = pdist2(Xi,X,dist);
    [~, idx] = sort(D);
    for j = 2:length(idx)
        if density(i) < density(idx(j))
            sigma(i) = pdist2(Xi,X(idx(j),:),dist);
            break;
        end
    end
    kp(i) = density(i)*sigma(i);
end
[~, kp_idx] = sort(kp,'descend');
kp_idx = kp_idx(1:keyPointsNum);
Xorg = X;
densOrg = density;
X = X(kp_idx,:);
density = densOrg(kp_idx);
%
%Repeat selection step for key points
while size(X,1)>clusterNum
    %
    %Key points
    sigma = inf(size(X,1),1);
    kp = inf(size(X,1),1);
    for i = 1:size(X,1)
        Xi = X(i,:);
        D = pdist2(Xi,X,dist);
        [~, idx] = sort(D);
        for j = 2:length(idx)
            if density(i) < density(idx(j))
                sigma(i) = pdist2(Xi,X(idx(j),:),dist);
                break;
            end
        end
        kp(i) = density(i)*sigma(i);
    end
    [~, new_idx] = sort(kp);
    D = pdist2(Xorg(kp_idx(new_idx(1)),:),Xorg(kp_idx,:));
    [~, I] = sort(D);
    % Use median value for next key point
    new_point = median([Xorg(kp_idx(new_idx(1)),:); Xorg(kp_idx(I(2)),:)]);
    I1 = find(kp_idx==kp_idx(new_idx(1)));
    I2 = find(kp_idx==kp_idx(I(2)));
    I = [I1; I2];
    kp_idx(I) = [];
    D = pdist2(new_point, Xorg);
    [~, I] = sort(D);
    kp_idx(end+1) = I(1);
    X = Xorg(kp_idx,:);
    density = densOrg(kp_idx);
    
end
kp_idx = kp_idx(:);

end

function kp_idx = keypointsComp2(X, dist, clusterNum, keyPointsNum)
% Description:
% Algorithm selects initial key points and removes points one-by-one until 
% clusterNum is reached. The point are removed based on relatively low density 
% and low density-based distance.  
%
% Inputs:
%            X - Input data set. Note: data must be complete.
%         dist - Selected distance metric. Default is: 'euc'.
%                Alternatives:
%                'sqe' - squared Euclidean distance
%                'euc' - Euclidean distance
%   clusterNum - The final number of clusters
% keyPointsNum - Initial number of key points
%
% Output:
%       kp_idx - Key point indexes
%
%Density
density = zeros(size(X,1),1);
for i = 1:size(X,1)
    Xi = X(i,:);
    D = pdist2(Xi,X,dist);
    [~, idx] = sort(D);
    knearests = X(idx(2:5),:);
    for j = 1:size(knearests,1)
        density(i) = density(i) + pdist2(Xi,knearests(j,:),dist);
    end
    density(i) = 1/density(i);
end
%
%Key points
sigma = inf(size(X,1),1);
kp = inf(size(X,1),1);
for i = 1:size(X,1)
    Xi = X(i,:);
    D = pdist2(Xi,X,dist);
    [~, idx] = sort(D);
    for j = 2:length(idx)
        if density(i) < density(idx(j))
            sigma(i) = pdist2(Xi,X(idx(j),:),dist);
            break;
        end
    end
    kp(i) = density(i)*sigma(i);
end
[~, kp_idx] = sort(kp,'descend');
kp_idx = kp_idx(1:keyPointsNum);
%
%Repeat selection step for key points
Xorg = X;
X = Xorg(kp_idx,:);
densOrg = density;
density = densOrg(kp_idx);
while size(X,1)>clusterNum
    %
    %Key points
    sigma = inf(size(X,1),1);
    kp = inf(size(X,1),1);
    for i = 1:size(X,1)
        Xi = X(i,:);
        D = pdist2(Xi,X,dist);
        [~, idx] = sort(D);
        for j = 2:length(idx)
            if density(i) < density(idx(j))
                sigma(i) = pdist2(Xi,X(idx(j),:),dist);
                break;
            end
        end
        kp(i) = density(i)*(sigma(i));
    end
    [~, new_idx] = sort(kp);
    I = find(kp_idx==kp_idx(new_idx(1)));
    kp_idx(I) = [];
    X = Xorg(kp_idx,:);
    density = densOrg(kp_idx);
    
end
kp_idx = kp_idx(:);

end

function params = parameters()
% Description: 
% Return initial parameters for test macro. 
%
% Output:
%     params - Parameters of test macro
%
params.datasets = {'S1'; 'S2'; 'S3'; 'S4'; 'S5D2'; 'S2D2'; 'O200'; 'O2000'; 'ecoli_2d'; 'iris_2d'; 'seeds_2d'};
params.correctResults = [15, 15, 15, 15, 5, 2, 5, 5, 3, 3, 3];
params.correctResults2 = [15, 15, 15, 15, 5, 2, 5, 5, 3, 2, 3];
params.probmiss = [0.00; 0.05; 0.10; 0.20];
params.maxClusters = 20;
params.treatMissingValidation = 'ads';
params.distance = 'euc';
params.neighborsNum = 2;
params.showProgression = true;
params.clusterInd = {@CalinskiHarabasz; @DaviesBouldin; @DaviesBouldin2; ...
    @GenDunn; @kCE; @PBM; @RayTuri; @Silhouette; @WB; @WemmertGancarski};

end


