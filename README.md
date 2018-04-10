# Entropy Clustering
This project addresses the problem of density-based clustering where clusters have variant densities and arbitrary shapes. Regarding the hardness of parameter regularization in existing methods, an entropy-based framework is given to generate clusters in two stages. Motivated by the fact that two clusters with ambiguous boundaries might be hardly discriminated, it presents the concept of density reachability to restrict the distance-based connectivity. By measuring entropy via k-NN distances in the first stage, a density order can be determined to glue the density-dominated cores into their adjacent dominators w.r.t. the density estimations. It highlights a heuristic algorithm in the second stage to effectively extract density peaks, such that the candidate peaks could enjoy a nice priority of low intra-cluster and high inter-cluster distances. The experimental study is conducted on dozens of data sets. The comparison results demonstrate that with simple tuning parameters, the proposed method shows substantial superiority to existing density-based clustering algorithms especially when the clusters exhibit various densities and arbitrary shapes.
## Invovled competitors
DBSCAN
DPC
KMEANS
## Results

Table 1: Comparison of the clustering accuracy.

Dataset|Metric| KMEANS| DBSCAN| PARMS |DPC| PARMS| EC| PARMS
--- | --- | --- | --- |--- |---| --- | --- | ---  
  -         |ACC |0.629 |0.977 |  -   |0.902 |  -    |0.992 |-
Compound    |AMI |0.556 |0.973 |1.456 |0.878 |0.368  |0.992 |3
  -         |ARI |0.683 |0.943 |/3    |0.844 |   -   |0.981 |-
  -         |ACC |0.802 |1.000 |  -   |1.000 |  -    |1.000 |-
Jain        |AMI |0.362 |1.000 |2.220 |1.000 |12.98  |1.000 |5
  -         |ARI |0.374 |1.000 |/14   |1.000 |   -   |1.000 |-
  -         |ACC |0.343 |1.000 |  -   |0.902 |  -    |1.000 |-
Spirals     |AMI |0.001 |1.000 |1.504 |0.878 |0.030  |1.000 |3
  -         |ARI |0.001 |1.000 |/3    |0.844 |   -   |1.000 |-
  -         |ACC |0.290 |0.900 |  -   |0.813 |  -    |0.983 |-
pathbased   |AMI |0.023 |0.766 |1.498 |0.597 |0.079  |0.953 |3
  -         |ARI |0.088 |0.702 |/6    |0.647 |   -   |0.931 |-
  -         |ACC |0.656 |0.710 |  -   |0.831 |  -    |0.819 |-
S3          |AMI |0.559 |0.481 |41884 |0.690 |3802   |0.673 |6
  -         |ARI |0.692 |0.678 |/108  |0.774 |   -   |0.755 |-
  -         |ACC |0.629 |0.600 |  -   |0.886 |  -    |0.910 |-
Seeds       |AMI |0.890 |0.458 |0.77  |0.663 |0.48   |0.697 |6
  -         |ARI |0.681 |0.350 |/10   |0.698 |   -   |0.749 |-
  -         |ACC |0.757 |0.877 |  -   |0.794 |  -    |0.891 |-
WDBC        |AMI |0.218 |0.478 |35    |0.293 |5.03   |0.507 |23
  -         |ARI |0.237 |0.656 |/10   |0.325 |   -   |0.611 |-
  -         |ACC |0.503 |0.305 |  -   |0.282 |  -    |0.690 |-
Segmentation|AMI |0.507 |0.356 |23.8  |0.156 |570    |0.612 |12
  -         |ARI |0.403 |0.280 |/26   |0.161 |   -   |0.553 |-
