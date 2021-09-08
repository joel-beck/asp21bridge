
<!-- README.md is generated from README.Rmd. Please edit that file -->

# asp21bridge

The goal of `asp21bridge` is to extend the
[lmls](https://gitlab.gwdg.de/asp21/lmls) package by implementing a
Markov Chain Monte Carlo Sampler with Ridge penalization.

## Underlying Model

The observation model in consideration is given by <!--  -->

``` math
\begin{aligned}
y_i \sim \mathcal{N} \left( \mathbf{x}_i^T \boldsymbol{\beta},\, \exp \left( \mathbf{ z}_i^T \boldsymbol{\gamma} \right)^2 \right) \qquad i = 1, \ldots n,
\end{aligned}
```

<!--  -->

where the *location* parameter *β* and the *scale* parameter *γ* are
themselves normally distributed with prior mean 0 and *hyperparameters*
*τ*<sup>2</sup> and *ξ*<sup>2</sup> for the prior variances.

## The Data

In the following examples the built-in simulated `toy_data` is used,
which consists of a column `y` representing a vector of observed values
and the explanatory variables `x1`, `x2`, `z1` and `z2`. The data is
constructed according to the mathematical model in the previous section,
where all explanatory variables predict the mean of `y` and only the
latter two model the variance.

``` r
library(asp21bridge)
head(toy_data)
## # A tibble: 6 x 5
##       y     x1    x2    z1    z2
##   <dbl>  <dbl> <dbl> <dbl> <dbl>
## 1 10.9  -0.207 0.194  5.41  2.62
## 2  6.91  1.28  1.42   4.53  3.10
## 3  9.40  2.08  0.891  5.07  4.64
## 4 10.5  -1.35  0.985  4.50  2.12
## 5  5.83  1.43  1.84   4.17  3.12
## 6  8.66  1.51  2.56   5.17  4.36
```

## Sampling Process

We first fit the frequentist `lmls` regression model and extend this
approach by adding the MCMC samples of the posterior distributions. For
later use, we also apply the MCMC sampler without Ridge - regularization
from the `lmls` package:

``` r
set.seed(1234)

fit <- toy_data %>%
  lmls(location = y ~ ., scale = ~ z1 + z2, light = FALSE) %>%
  mcmc(nsim = 1000) %>%
  mcmc_ridge(num_sim = 1000)
```

## Numerical Analysis

The `asp21bridge` package contains various tools to analyze the sampling
results both numerically and graphically.

A first quick overview can be gained by the generic `summary()` function
with specification of the `type` argument:

``` r
summary(fit, type = "mcmc_ridge")
## 
## Call:
## lmls(location = y ~ ., scale = ~z1 + z2, data = ., light = FALSE)
## 
## Pearson residuals:
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## -3.29800 -0.38660  0.11770 -0.01354  0.57410  2.54600 
## 
## Location coefficients (identity link function):
##             Mean      2.5%       50%  97.5%
## beta_0 -0.003064 -0.152295 -0.001275  0.132
## beta_1 -2.001501 -2.007618 -2.001631 -1.995
## beta_2 -1.004306 -1.017451 -1.004198 -0.991
## beta_3  1.001918  0.984266  1.001849  1.022
## beta_4  2.008085  2.002211  2.007997  2.014
## 
## Scale coefficients (log link function):
##            Mean    2.5%     50%  97.5%
## gamma_0  0.9257 -0.1398  1.0384  1.982
## gamma_1 -1.1528 -1.3285 -1.1598 -0.956
## gamma_2  0.9403  0.7134  0.9438  1.185
## 
## Residual degrees of freedom: 42 
## Log-likelihood: 32.28 
## AIC: -48.57 
## BIC: -33.27
```

A more comprehensive list of the estimated parameter values is provided
by the `summary_complete()` function:

``` r
summary_complete(fit)
## # A tibble: 10 x 6
##    Parameter `5% Quantile` `Posterior Mean` `Posterior Median` `95% Quantile`
##    <chr>             <dbl>            <dbl>              <dbl>          <dbl>
##  1 beta_0           -0.123         -0.00306           -0.00127          0.107
##  2 beta_1           -2.01          -2.00              -2.00            -2.00 
##  3 beta_2           -1.01          -1.00              -1.00            -0.994
##  4 beta_3            0.987          1.00               1.00             1.02 
##  5 beta_4            2.00           2.01               2.01             2.01 
##  6 gamma_0          -0.106          0.926              1.04             1.84 
##  7 gamma_1          -1.31          -1.15              -1.16            -0.984
##  8 gamma_2           0.750          0.940              0.944            1.13 
##  9 tau^2             3.14           4.01               3.94             5.09 
## 10 xi^2             15.9           51.0               37.4            132.   
## # ... with 1 more variable: Standard Deviation <dbl>
```

Since the results are embedded in a `data frame`, the usual methods of
data frame manipulation allow for a convenient analysis even for high
dimensional parameter vectors. A particularly interesting use case is
the comparison of the coefficient estimates from the following three
approaches:

-   Point estimates based on Maximum Likelihood Estimation (`mle`
    column)

-   Posterior Mean estimates from MCMC sampling without regularization
    (`mcmc` column)

-   Posterior Mean estimates with Ridge penalty (`mcmc_ridge` column)

Both the `mcmc()` function as well as the `mcmc_ridge()` function simply
add their results to the existing model. Therefore all relevant
information is contained in the `fit` object from above:

``` r
summary_complete(fit) %>%
  dplyr::transmute(Parameter, mcmc_ridge = `Posterior Mean`) %>%
  dplyr::filter(stringr::str_detect(Parameter, pattern = "beta|gamma")) %>%
  dplyr::mutate(
    mcmc = summary_complete(fit$mcmc)$`Posterior Mean`,
    mle = c(coef(fit)$location, coef(fit)$scale),
    truth = c(0, -2, -1, 1, 2, 0, -1, 1)
  )
## # A tibble: 8 x 5
##   Parameter mcmc_ridge     mcmc     mle truth
##   <chr>          <dbl>    <dbl>   <dbl> <dbl>
## 1 beta_0      -0.00306 -0.00762 -0.0217     0
## 2 beta_1      -2.00    -2.00    -2.00      -2
## 3 beta_2      -1.00    -1.00    -1.00      -1
## 4 beta_3       1.00     1.00     1.00       1
## 5 beta_4       2.01     2.01     2.01       2
## 6 gamma_0      0.926    1.17     1.44       0
## 7 gamma_1     -1.15    -1.18    -1.31      -1
## 8 gamma_2      0.940    0.920    0.999      1
```

The estimates for *β*<sub>1</sub> up to *β*<sub>4</sub> are identical in
the three models and almost coincide with the true data generating
values. This behaviour is somewhat expected, since the *β* vector is
drawn from a closed form multivariate normal distribution with
independent samples across iterations of the MCMC sampler.

The results for the scale parameter **γ** show a greater variation
between the models and are therefore more interesting to analyze. The
Maximum Likelihood approach is able to recover the true / data
generating value of *γ*<sub>2</sub>, while the MCMC sampler with Ridge
penalty performs best for estimating *γ*<sub>1</sub>. All three models
fail to identify the correct intercept parameter *γ*<sub>0</sub>, which,
however, is often of minor interest.

Simulations across different data sets and sample sizes show similar
patterns: If the data is generated according to the mathematical model
introduced in the first section, i.e. the estimation model using Maximum
Likelihood is correctly specified, the `mle` and the `mcmc_ridge`
estimates tend to be closest to the true values. If, however, the
underlying model used for simulation deviates from the Gaussian
location-scale regression model, the `mle` results are, unsurprisingly,
worse than both of the more flexible MCMC methods.

The Markov Chain Monte Carlo samples for **γ** based on the Metropolis -
Hastings algorithm implemented in the `mcmc_ridge()` function are
strongly correlated, even though the acceptance rate of the proposed
values is reasonable:

``` r
fit$mcmc_ridge$acceptance_rate
## [1] 0.408
```

This value of roughly 41% indicates that the variance of the proposal
distribution is chosen appropriately.

## Graphical Analysis

The building blocks for monitoring the convergence of the posterior
chains as well as the autocorrelations are the functions
`diagnostic_plots()` for a single Markov Chain and `mult_plot()` for
combining multiple chains.

A quick overview can be gained by collecting the corresponding trace
plots for all posterior coefficients.

``` r
mult_plot(fit, type = "trace", free_scale = TRUE, latex = TRUE)
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />

The trace plots for the **β** coefficients indicate convergence and
confirm the stable posterior estimates from the previous section.

The prior variance parameters *τ*<sup>2</sup> and *ξ*<sup>2</sup> are
strictly positive and thus more adequately displayed on the logarithmic
scale:

``` r
samples <- fit$mcmc_ridge$sampling_matrices
mult_plot(
  samples = list(samples$location_prior, samples$scale_prior),
  type = "both", log = TRUE, latex = TRUE
)
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

Just like the **β** samples, these are drawn from closed form full
conditional distributions.

Finally, we focus on the samples for the scale parameter *γ*. Here, the
trace plots do not indicate convergence to the posterior distribution.
Thus, we increase the number of simulations to 10000, which also allows
the inclusion of a Burn - In Phase as well as a thinning factor while
still maintaining a sufficient sample size to estimate posterior
quantities:

``` r
set.seed(4321)
fit <- fit %>%
  mcmc_ridge(num_sim = 10000)
```

We drop the first 1000 samples for *γ*<sub>2</sub> and show the three
most common diagnostic plots for the thinned sample:

``` r
samples <- fit$mcmc_ridge$sampling_matrices

samples$scale[, 3, drop = FALSE] %>%
  burnin(num_burn = 1000) %>%
  thinning(freq = 30) %>%
  diagnostic_plots(lag_max = 30, latex = TRUE)
```

<img src="man/figures/README-unnamed-chunk-11-1.png" width="100%" />

Even a thinning factor of 30, i.e. only every 30th observation is kept
in the sample, does not get rid of the autocorrelation. Yet, the
posterior density seems to be approximately normal and the trace plot
indicates a decent exploration of the posterior space.

These findings serve as a great starting point for a more in-depth
analysis, which is beyond the scope of this brief introduction to the
`asp21bridge` package. Full descriptions of all functions, their
arguments and common applications can be found in the documentation.
