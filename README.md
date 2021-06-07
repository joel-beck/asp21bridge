
<!-- README.md is generated from README.Rmd. Please edit that file -->

# asp21bridge

The goal of `asp21bridge` is to extend the `lslm` package by implementing
a Markov Chain Monte Carlo Sampler with Ridge penalization.

## Underlying Model

The observation model in consideration is given by
<!--  -->

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
and the explanatory variables `x1`, `x2`, `z1` and `z2`. Here, all
explanatory variables predict the mean of `y` and only the latter two
model the variance.

``` r
library(asp21bridge)
head(toy_data)
## # A tibble: 6 x 5
##       y     x1    x2    z1    z2
##   <dbl>  <dbl> <dbl> <dbl> <dbl>
## 1  29.0 -0.207 0.194  5.41  2.62
## 2  23.3  1.28  1.42   4.53  3.10
## 3  30.8  2.08  0.891  5.07  4.64
## 4  18.7 -1.35  0.985  4.50  2.12
## 5  19.8  1.43  1.84   4.17  3.12
## 6  23.7  1.51  2.56   5.17  4.36
```

## Sampling Process

We first apply the usual regression model based on Maximum - Likelihood
estimation implemented in the `lslm` model and enhance this approach by
adding the MCMC samples of the posterior distributions:

``` r
set.seed(1234)

fit <- toy_data %>%
  lslm(location = y ~ ., scale = ~ z1 + z2, light = FALSE) %>%
  gibbs_sampler(num_sim = 1000)
```

## Numerical Analysis

The `asp21bridge` package contains various tools to analyze the sampling
results both numerically and graphically.

A first quick overview can be gained by the generic `summary()`
function:

``` r
summary(fit, type = "mcmc_ridge")
##
## Call:
## lslm(location = y ~ ., scale = ~z1 + z2, data = ., light = FALSE)
##
## Pearson residuals:
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
## -3.07000 -0.55750  0.06667  0.00000  0.59940  2.27500
##
## Location coefficients (identity link function):
##           Mean    2.5%     50%  97.5%
## beta_0  0.5468 -0.9426  0.5800  1.971
## beta_1  2.2240  1.8702  2.2243  2.543
## beta_2 -4.0308 -4.2951 -4.0325 -3.739
## beta_3  4.9699  4.6681  4.9763  5.252
## beta_4  0.9182  0.6892  0.9198  1.149
##
## Scale coefficients (log link function):
##              Mean      2.5%       50% 97.5%
## gamma_0  0.104034 -1.170869 -0.003655 1.550
## gamma_1  0.004717 -0.228834  0.006573 0.168
## gamma_2 -0.024053 -0.208590 -0.035959 0.189
##
## Residual degrees of freedom: 42
## Log-likelihood: -68.83
## AIC: 153.7
## BIC: 169
```

A more comprehensive list of the estimated parameter values is provided
by the `summary_complete()` function:

``` r
results <- summary_complete(fit)
results
## # A tibble: 10 x 6
##    Parameter `5% Quantile` `Posterior Mean` `Posterior Median` `95% Quantile`
##    <chr>             <dbl>            <dbl>              <dbl>          <dbl>
##  1 beta_0         -0.711            0.547              0.580           1.72
##  2 beta_1          1.95             2.22               2.22            2.50
##  3 beta_2         -4.26            -4.03              -4.03           -3.79
##  4 beta_3          4.70             4.97               4.98            5.21
##  5 beta_4          0.724            0.918              0.920           1.11
##  6 gamma_0        -1.03             0.104             -0.00366         1.49
##  7 gamma_1        -0.208            0.00472            0.00657         0.158
##  8 gamma_2        -0.192           -0.0241            -0.0360          0.159
##  9 tau^2           0.00530          0.0153             0.0117          0.0309
## 10 xi^2            0.0545           0.212              0.140           0.552
## # ... with 1 more variable: Standard Deviation <dbl>
```

Since the results are embedded in a `data frame`, the usual methods of
data frame manipulation allow for a convenient analysis even for high
dimensional parameter vectors. For illustration purposes we can easily
sort all scale parameters by their standard deviation:

``` r
results %>%
  dplyr::select(Parameter, `Standard Deviation`) %>%
  dplyr::filter(stringr::str_detect(Parameter, pattern = "gamma")) %>%
  dplyr::arrange(dplyr::desc(`Standard Deviation`))
## # A tibble: 3 x 2
##   Parameter `Standard Deviation`
##   <chr>                    <dbl>
## 1 gamma_0                  0.764
## 2 gamma_2                  0.115
## 3 gamma_1                  0.106
```

## Graphical Analysis

The building blocks for monitoring the convergence of the posterior
chains as well as the autocorrelations are the functions
`diagnostic_plots()` for a single Markov Chain and `mult_plot()` for
combining multiple chains.

First, we analyze the convergence of the location coefficients. Note
that the posterior samples are saved in the list entry
`mcmc_ridge$sampling_matrices`:

``` r
samples <- fit$mcmc_ridge$sampling_matrices

mult_plot(
  samples = samples$location, type = "both",
  free_scale = TRUE, latex = TRUE
)
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

We can observe stable time plots indicating convergence and a sufficient
exploration of the posterior space as well as approximate normal
posterior distributions for *β*<sub>1</sub> up to *β*<sub>4</sub>.

While time plots, density plots and autocorrelation plots can always be
displayed separately, it is often convenient to combine all of them in a
single plot. This can be achieved with the `diagnostic_plots()` function
as demonstrated below. Since the variance parameter *ξ*<sup>2</sup> is
strictly positive, we choose a logarithmic scale:

``` r
diagnostic_plots(samples = samples$scale_prior, log = TRUE, latex = TRUE)
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />

All three plots confirm the validity of the posterior distribution
estimates. In this case there is neither a thinning procedure nor a
Burn-In phase necessary, which are implemented by the functions
`thinning()` and `burnin()`.
