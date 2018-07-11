## **Bayesian integration w/ tmbstan**
In this example you learn:

* How to use tmbstan to do MCMC sampling of a TMB model using the software
  tmbstan. 
* Mix integration the integration techniques of the Laplace approximation
  and Markov chain Monte Carlo (MCMC) integration on the same model with no
  additional effort.
* Test the accuracy of the Laplace approximation to the marginal likelihood

***
#### Description of the practical situation
Here the model itself is not as important, but it is a binomial GLMM with
three levels of crossed random effects. This particular analysis comes from
the paper Monnahan and Kristensen (2018), where you can find a more
detailed description of the model and analysis, as well as further
resources.

Here we simply note this model is setup to be Bayesian, meaning it has
explicit priors and Jacobian adjustments. Further, the three hypervariances
are in log space (yearInterceptSD, plantInterceptSD, plantSlopeSD).

Instead, we are interested in testing the the Laplace approximation (LA)
used by TMB which approximates the marginal likelihood (i.e., the
likelihood with the random effects integrated out). Typically, there is no
easy way to test how valid the approximation is. Here, we demonstrate how
to do this by integrating the model in two different ways: (1) with MCMC
sampling for all parameters, and (2) with MCMC for fixed effects but LA for
random effects. That is, the MCMC algorithm samples from the approximated
marginal posterior. If the LA is accurate, then there will be no difference
in the posterior distribution for the *fixed effects*. We cannot directly
compare the random effects.

Thus, if the distribution is the same we have confidence that the LA is
accurate. If not, this suggests an inaccuracy in the LA and caution should
be taken when using it for inference. In the case of this model, we quickly
see that the LA version produces very different results for the fixed
effects.

Another interesting note is that while the LA version mixes better in the
sense that it produces higher effective samples, it takes much longer to
run and thus the full MCMC integration should be prefered from an
efficiency standpoint as well (see table S2 in Monnahan and Kristensen
2018).


***
#### References
Monnahan, C. C. and K. Kristensen. 2018. No-U-turn sampling for fast
Bayesian inference in ADMB and TMB: Introducing the adnuts and tmbstan R
packages. Plos One 13:e0197954.
