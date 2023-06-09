# Randomized Feature Selection Based Semi-Supervised Latent Dirichlet Allocation for Microbiome Analysis


R code for  Randomized Feature Selection Based Semi-Supervised Latent Dirichlet Allocation for Microbiome Analysis
## Abstract:
Health and disease are fundamentally influenced by microbial communities and their genes (the microbiome).
An in-depth analysis of microbiome structure that enables the classification of individuals based on their health can be crucial in enhancing diagnostics and treatment strategies to improve the overall well-being of an individual.

In this paper, we present a novel semi-supervised methodology known as Randomized Feature Selection based Latent Dirichlet Allocation (RFSLDA) to study the impact of the gut microbiome on a subject's health status. Since the data in our study consists of fuzzy health labels, which are self-reported, traditional supervised learning approaches may not be suitable. As a first step, based on the similarity between documents in text analysis and gut-microbiome data, we employ Latent Dirichlet Allocation (LDA), a topic modeling approach which uses microbiome counts as features to group subjects into relatively homogeneous clusters, without invoking any knowledge of observed health status (labels) of subjects. We then leverage information from the observed health status of subjects to  associate these clusters with the most similar health status making it a semi-supervised approach. Finally, a feature selection technique is incorporated into the model to improve the overall classification performance.

### Results:
The proposed method provides a semi-supervised topic modelling approach that can help handle the high dimensionality of the microbiome data in association studies. Our experiments reveal that our semi-supervised classification algorithm is  effective and efficient in terms of high classification accuracy compared to popular supervised learning approaches like SVM and multinomial logistic model.

### Conclusion:
The RFSLDA framework is attractive because it (i) enhances clustering accuracy by identifying key bacteria types as indicators of health status, (ii) identifies key bacteria types \textit{within each group} based on  estimates of the proportion of bacteria types within the groups,  and (iii) computes a measure of within-group similarity to identify highly similar subjects in terms of their health status. 


## Dependencies:
Code was developed using R version 4.2.3.
## Description: 
1) **MicrobiomeData_109.rda**: Complete data .
2) **Data.rda**: Feature data( Bacteria counts) .
3) **actual.rda**: Observed health labels .
4) **RFSLDA.R**: performs RFSLDA on Data  .
5) **lda_final_iters_fromtop50.RData**: RFSLDA results saved after running the RFSLDA for R=50 iterations. 
6) **RFSLDAValidation.R**: Compares the RFSLDA method with popular supervised learning methods like Support vector machine and multinomial logistic model using a train test validation.
RFSLDA.pdf: Preprint of the paper (Under submission):  All the results obtained here are presented in the paper.






