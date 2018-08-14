# Bayesian analysis of individual treatment response (ITR) curves on EHR time series

The main goal of this GSoC project is to speed up the current ITR [1] implementation, which was purely written in R, by using Rcpp to code the core inference part. ITR is a Bayesian approach that posits hierarchical priors (i.e. DPM) on the parameters of the treatment response curves, which results at Bayesian estimators at both group level and individual level. The reason of using DPM prior is because it allows individual parameterization but also can borrow strength from others to estimate individual parameters. Here, the treatment response model assumes a Gaussian process regression model for the baseline progression of a patient outcome variable (i.e. with no exposure to treatments), and assumes the response to multiple treatments are additive to the baseline progression. The response to a treatment (i.e. treatment effect) is modeled in a U-shape curve with five free parameters specifying the rate that the effect climbs up at the beginning, when and where the effect reaches its peak, the rate that the effect goes down afterwards and the final value the effect stays at permanently. 
The main MCMC inference is done by Gibbs sampling since most of the variables are unconstrained and their priors are assumed to be conjugate. Constrained variables, such as the climbing-up rate and going-down rate, are first transformed to be unconstrained and use Jacobian adjustment to posit priors; for these non-conjugated cases, Metropolis-Hastings samplers are developed. The sampling inference for DPM is done using the truncated stick-breaking process that was developed by [2]. 

In this package, inference of ITR can be implemented using three models: DPM-ITR estimating treatment responses at both group level and individual level; DP-ITR estimating treatment responses only at group level; ITR estimating treatment responses at individual level. A [demo.R] is included for demonstating the three implementations on a simulated data.

What have done.
- Rcpp implementations of the three MCMC inference method: DPM-ITR, DP-ITR, ITR
- A R function is included in the package for simulating longitudinal data that include
  - irregularly sampled and sequential exposures to two different treatment, one of which has positive effects on the outcome and the other has negative effects.
  - irregularly sampled observations of the outcome variable.
  - heterogenous baseline progression and treatment response.
- A R function is also provided for plotting the estimated baseline progression and treatment response curves.

What to do.
- A sample real data set will be included and another demo will be provided for demonstrating how the package can be configured in real experiments.
- Inference on newly test data will be integrated and various evaluation metrics will also be provided.
- A R plot function will be included for visulizations like Figure 15 in [1].

[1]  Xu, Y., Xu, Y. and Saria, S., 2016. A Bayesian Nonparametic Approach for Estimating Individualized Treatment-Response Curves. arXiv preprint arXiv:1608.05182. [accepted as minor revision by Journal of Machine Learning Research]
[2] Ishwaran, H. and James, L. F. (2001). Gibbs sampling methods for stick-breaking priors. Journal of the American Statistical Association 96. 


