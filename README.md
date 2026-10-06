# CoCluster

**CoCluster** is an R package designed to simultaneously cluster (co-cluster) the rows and columns of data matrices with native support for ordinal data.

Standard clustering algorithms group only observation rows. In high-dimensional questionnaire data—such as quality-of-life surveys where responses are given on ordinal Likert scales—co-clustering is far more powerful. It reduces the entire dataset into a concise summary of "types of respondents" and "types of questions" by finding homogeneous, interacting blocks of data.

This package leverages the **Latent Block Model (LBM)** equipped natively with the **Combination of Uniform and shifted Binomial (CUB)** observation model, tailored for ordinal data to account for both deliberate human choices and random uncertainty, and can easily be extended to use other models. It estimates the optimal co-clustering using the **SEM-Gibbs** algorithm and can automatically search for the optimal number of clusters using a greedy surface search algorithm and the **Integrated Completed Likelihood (ICL)** criterion.

## To dos

un package R (en format .tar.gz installer sur Mac)