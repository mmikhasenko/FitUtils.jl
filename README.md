# FitUtils.jl

The package is a heavy extension of the [AlgebraPDF.jl](https://github.com/mmikhasenko/AlgebraPDF.jl) that provide convenient interface for fitting with `Optim.jl` and `Minuit` (via the `PyCall` interface).

## How to Fit
A set of data can be fit by minimizing negative log-likelihood function.
A simple method `fit` combines three steps:
1. Create the negative log-likelihood function
2. Call minimization algorithm (gradient decent and error evaluation)
3. Build a summary combining parameters with the errors, and the model, updated to new values (best model)

```julia
d = Normalized(FGauss((μ=0.6, σ=1.5)), (-5, 5))
data = randn(1000)
fit_summary = fit(d, data)

fit_summary.parameters # (μ = -0.023040771466510682, σ = 0.9847923919821122)
fit_summary.measurements # (μ = -0.023 ± 0.031, σ = 0.985 ± 0.022)
# 
plot(fit_summary.best_model) # plot just best model
plot(data, fit_summary.best_model, bins=30, lw=3, lab="model") # plot the fitted model together with the data
```

The negative log-likelihood function is a subtype of `AbstactFunctionWithParameters`.
An optional argument is a value which is used instead of logarithm in case the `model(x) < 0` (might happen during the minimization process).
```julia
nll = NegativeLogLikelihood(model, data, nagativepenatly=-1e4)
nll = minussum(log(model))
```
The second method does the same, but perhaps more intuitive.

There are two options to pass to the fitting faction for minimization algorithm.
```julia
fs1 = fit(model, data, MigradAndHesse(errordef=1/2)) # calls Minuit from iminuit 
fs2 = fit(model, data, BFGSApproxHesse()) # calls BFGS of Optim package
``` 
The first method is found to be extremely stable and fast (despite finite-diff derivatives),
therefore, used by default. The minimization algorithm is a predecessor of `BFGS` but equipped
by objective stopping criterion [see old Minuit paper](https://www.sciencedirect.com/science/article/abs/pii/0010465575900399).
The second method calls the `BFGS()` algorithm implemented in `Optim` that has a great convergence,
however, noted to fail often doing the line search (HagerZhang) with the error `isfinite(phi_c) && isfinite(dphi_c)`. It happens particularly when parameters run away of the reasonable range and numerical value of
the normalization function becomes 0.
The errors are computed after minimization step using either `HESSE` algorithm of Minuit,
or approximated hessian matrix attached to the state of Optim.

The summary of the fit result returned by the fit function is a named tuple containing
 * `parameters` - the values of parameters that minimize the NLL,
 * `measurements` - the values of parameters with errors
 * `best_model` - the input model with parameters updated
 * `nll` - the value of NLL in the minimum, and
 * `fit_result` - the object returned by the minimizer. Can be used, e.g. to get the covariance matrix.

For the model constructed as a sum of normalized functions, an extended maximum likelihood fit (ENLL includes a Poisson term for normalization) can be performed.
Here is an example:
```julia
signal = Normalized(abs2(FBreitWigner((m=0.77,G=0.15))), (0.28, 1.1))
backgr = Normalized(FExp((τ=-1.3,)), (0.28, 1.1))

model = signal + backgr # calls +([signal, backgr], (α1=1.0, α2=1.0))

data = generate(model, 1000)

let binning = range(lims(model)..., length=80)
    stephist(data, bins=binning )
    plot!(model, scaletobinneddata(1000/(1+1),binning)) # 1+1 is the current normalization of the model
end

model_init = updatepars(model, (α1=400, α2=600))
enll = Extended(NegativeLogLikelihood(model_init, data))
# 
fit_summary = fit(enll)
fit_summary.measurements
# returns
# (m = 0.7777 ± 0.0063,
#  G = 0.17 ± 0.026,
#  τ = -1.31 ± 0.36,
#  α1 = 540.0 ± 67.0,
#  α2 = 460.0 ± 66.0)

let binning = range(lims(model)..., length=80)
    stephist(data, bins=binning)
    plot!(fit_summary.best_model, scaletobinneddata(binning)) # 1+1 is the current normalization of the model
end
```
