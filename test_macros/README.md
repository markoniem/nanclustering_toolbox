#### Experiments related to the toolbox ####

- **part0_toolboxdemo:** demonstrates the cluster validation pipeline starting from the pre-processing phase and ending to the validation phase.
- **part0_extended_pdist2:** uses the expected squared Euclidean distance for computing distances between observations and centroids.
- **part0_useParallel:** Performs weighted K-spatial-medians clustering using parallelled computation. 
- **part0_useGPU:** Performs weighted K-spatial-medians clustering using GPU computation. 
- **partA_evaluate_dist_functions:** compares distance computation functions (ADS, PDS, ESD, EED, ICkNNI (k=5), kNNI (k=5), IST_MC, ICkNNI_EXP (k=5), and kNNI_EXP (k=5)) using data with missing values. Results are available in [**partA_Results.xlsx**](partA_Results.xlsx) 
- **partB_evaluate_cluster_indices:** compares ten internal cluster validation indices using K-spatial-median clustering with ADS and K-spatial-median clustering with EED-ADS. Results are in [**partB_Results1**](partB_Results1/) and in [**partB_Results2**](partB_Results2/).
- **partB_estimation_in_indices:** benefits distance estimation in actual cluster validation indices. Clustering is performed using K-spatial-median with ADS and K-spatial-median with EED-ADS. Results are in [**partB_Results1b.xlsx**](partB_Results1b.xlsx) and in [**partB_Results2b.xlsx**](partB_Results2b.xlsx).  
- **partB_determine_quality_of_clustering:** uses three external validation indices (ACC, ARI, and NMI) to measure performances of K-spatial-median clustering with ADS and K-spatial-median clustering with EED-ADS. Results are in [**partB_Results3.xlsx**](partB_Results3.xlsx).        
- **partB_key_point_selection:** uses key point selection algorithm in initialization of K-spatial-medians clustering. Cluster validation results are in [**partB_key_point_results**](partB_key_point_results/).
- **partB_rmse_results:** computes RMSE errors between real centroids used in initialization of K-spatial-medians clustering and centroids obtained with K-means++ initialization of K-spatial-medians clustering. Experiments are performed using clustering based on ADS, EED, and EED-ADS distance computation methods. Results are available in [**partB_rmse_results**](partB_rmse_results/). See detailed description of experiments in [**Distance_estimation_in_clustering.pdf**](Distance_estimation_in_clustering.pdf).       
- **partC_validate_multiDim_data:** utilizes cluster validation for multidimensional data with missing values. Experiments are performed for data sets with various dimensionalities (M=10-100) and cluster overlap (dc=0.6-0.9). Results are in [**partC_Results.xlsx**](partC_Results.xlsx).  


