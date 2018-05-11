# SPICE
 Stochastic Parameter Inference with the Cross-Entropy Method: 


Examples contained within /models/:
* Lotka-Volterra Model
* Yeast Polarization Model
* Schl&ouml;gl system
* Toggle Switch (with modified src files to account for model specific behaviour)

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

# ---- Bounds on the initial (log) search space for parameters
bounds = [[1e-6,1e-6,1e-6],
		[10.0,10.0,10.0]]


# ---- We specify seperable mass-action propensity functions (acting on a path p, containing as attributes the species x, parameters θ, propensity vector a, time t)
function F(p)
    H(p)
    p.a .*= p.θ
end

function H(p)
    p.a[1] = p.x[1]
    p.a[2] = p.x[1] * p.x[2]
    p.a[3] = p.x[2]
end

# ---- Specify the model to SPICE, passing information that the observables in the data are called X1 and X2.
model = Model(x0, F, H, ν, bounds, obsname=[:X1,:X2])

# ---- We can tweak the algorithm parameters below
# e.g., ssa = :Tau (or :Direct)
cem = CEM(ssa=:Tau, nElite = 10, nRepeat=1, nSamples=1000, maxIter=250, mSamples=20000, shoot=false, splitting=false, sampling=:log, tauOpt=TauOpt(ϵ=0.1))

# ---- Initialise the system, and point it to the folder containing the data
system = System(model, "./data/5", routine=cem)

# ---- run 1 run of the algorithm, appending output results with the string "lotka"
estimate(system, 1, "lotka")
```

All outputs are sent to the "/results" folder. Given M parameters during estimation, the output CSV files contain M columns with the final (log) parameter estimates, M with their associated variances, and the final column is the CPU time. Further iteration-by-iteration diagnostics are provided in "/results/traces".
