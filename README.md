# SPICE
 Stochastic Parameter Inference with the Cross-Entropy Method:

Examples contained within /models/:
* Lotka-Volterra Model
* Yeast Polarization Model
* Schl&ouml;gl system
* Toggle Switch (with modified src files)

To run a model from the command line:

`$ julia model.jl` 

where model.jl is replaced with the name of the model of interest.

## Requirements:
* Julia (v0.6+) [see here](https://github.com/JuliaLang)

## Packages required:
* Distributions
* Sobol
* DataFrames
* Compat
* StatsBase

To add packages within Julia, use -
```julia
Pkg.add("Name") 
```

### Advanced Example
Here we consider the Lotka-Volterra predator-prey model as an example:

		X1 -> X1 + X1
		X1 + X2 -> X2 + X2
		X2 -> 0

We can specify this mass-action model within Julia in the following way:

```julia
# include the SPICE module
include("../../src/SPICE.jl")
using SPICE

# ---- add the State-change matrix
ν = [1 0; -1 1; 0 -1]

# ---- add the Initial state for the species
x0 = [50,50]

# ---- Bounds on the initial search space for parameters
bounds = [[0.0,0.0,0.0], 
		[1.0,1.0,1.0]]

# ---- We specify separable mass-action propensity functions (acting on species x, parameters θ, propensity vector a, time t)
function F(x, θ, a, t)
    H(x,θ,a,t)
    a .*= θ
end

function H(x, θ, a, t)
    a[1] = x[1]
    a[2] = x[1] * x[2]
    a[3] = x[2]
end

# ---- Specify the model to SPICE, passing information that the observables in the data are called X1 and X2.
model = Model(x0, F, H, ν, bounds, obsname=[:X1,:X2])

# ---- We can tweak the algorithm parameters below
# e.g., ssa = :Tau (or :Direct)
cem = CEM(ssa=:Tau, nElite = 10, nRepeat=1, nSamples=1000, maxIter=250, mSamples=20000, shoot=true, splitting=false, tauOpt=TauOpt(ϵ=0.1))

# ---- Initialise the system, and point it to the folder containing the data
system = System(model, "./data/5", routine=cem)

# ---- run 1 run of the algorithm, appending output results with the string "lotka"
estimate(system, 1, "lotka")
```
