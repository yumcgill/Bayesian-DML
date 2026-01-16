
This repository contains the R code used to implement the simulation studies in the paper:

[**A Scalable Bayesian Double Machine Learning Framework for High-Dimensional Causal Estimation, with Application to Racial Disproportionality Assessment**](https://arxiv.org/abs/arXiv:2502.07695)  by Yu Luo, Vanessa McNealis, and Yijing Li

The code reproduces the simulation results for both binary and continuous treatment settings using the proposed Bayesian double machine learning framework.

## Repository Structure
The repository is organized into two main folders corresponding to the two simulation examples in the paper.

## Example 1: Binary treatment
The ```Example 1/``` folder contains code for simulations with a binary treatment variable.

### Files
+ ```Neyman_EL_Ex1.R```
 Implements MCMC using Algorithm 2 with $K=2$ based on the Bayesian Empirical Likelihood (EL) method.

+ ```Neyman_ETEL_Ex1.R```
 Implements MCMC using Algorithm 2 with $K=2$ based on the Bayesian Exponentially Tilted Empirical Likelihood (ETEL) method.

+ ```Neyman_HD_Ex1.R```
 Implements MCMC using Algorithm 2 with $K=2$ based on the Bayesian High-Dimensional (HD) method.

+ ```runEx1.R```
 Runs the simulation study for:
  + Sample sizes: $n=50$ and $n=200$
  + Number of replications: 1,000

## Example 2: Continuous Treatment
The ```Example 2/``` folder contains code for simulations with a continuous treatment variable.

### Files 
+ ```Neyman_EL_Ex2.R```
Implements MCMC using Algorithm 2 with $K=2$ based on the Bayesian Empirical Likelihood (EL) method.

+ ```Neyman_ETEL_Ex2.R```
Implements MCMC using Algorithm 2 with $K=2$ based on the Bayesian Exponentially Tilted Empirical Likelihood (ETEL) method.

+ ```Neyman_HD_Ex2.R```
Implements MCMC using Algorithm 2 with $K=2$ based on the Bayesian High-Dimensional (HD) method.

+ ```runEx2.R```
 Runs the simulation study for:
  + Sample size: $n=40$
  + Number of replications: 1,000


## Note

+ All methods are implemented using MCMC according to Algorithm 2 described in the paper.
+ The simulation scripts, ```runEx1.R``` and ```runEx2.R```, source the corresponding method files and generate summary results.

