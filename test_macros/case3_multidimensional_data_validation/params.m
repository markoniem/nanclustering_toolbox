function params = params()
% Description: 
% Return initial parameters for test macros. See more detailed
% instructions of parameters from original macros.  
%
% Output:
%     params - Parameters of test macro
%
params.datasets = {'K15_Nk400_M10_SPEHERES_dc0.9.mat'; 'K15_Nk400_M10_SPEHERES_dc0.8.mat';
                   'K15_Nk400_M10_SPEHERES_dc0.7.mat'; 'K15_Nk400_M10_SPEHERES_dc0.6.mat'
                   'K15_Nk400_M50_SPEHERES_dc0.9.mat'; 'K15_Nk400_M50_SPEHERES_dc0.8.mat';
                   'K15_Nk400_M50_SPEHERES_dc0.7.mat'; 'K15_Nk400_M50_SPEHERES_dc0.6.mat';
                   'K15_Nk400_M100_SPEHERES_dc0.9.mat'; 'K15_Nk400_M100_SPEHERES_dc0.8.mat';
                   'K15_Nk400_M100_SPEHERES_dc0.7.mat'; 'K15_Nk400_M100_SPEHERES_dc0.6.mat'};             
params.probmiss = [0.00; 0.05; 0.10; 0.20];
params.maxClusters = 20;
params.clustMethod = @kcentroids;
params.distance = 'euc';
params.replicates = 100;
params.useprevCent = true;
params.initcrit = 'kmeans++';
params.pipelined = true;
params.showProgression = true;
params.distFuncInd = 'ads';
params.clusterInd = {@CalinskiHarabasz; @DaviesBouldin; @DaviesBouldin2; ... 
   @GenDunn; @kCE; @PBM; @RayTuri; @Silhouette; @WB; @WemmertGancarski};



