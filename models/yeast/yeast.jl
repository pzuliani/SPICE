include("../../src/SPICE.jl")
using SPICE

# ---- State-change matrix
ν = [[1 0 0 0 0 0 0];[0 0 0 0 0 0 0];[0 1 1 0 0 0 0];[1 0 0 0 0 0 0];[0 0 0 0 1 1 0];[0 0 0 0 0 0 1];[0 0 0 1 0 0 0];[0 0 1 0 0 0 0]]-
[[0 0 0 0 0 0 0];[1 0 0 0 0 0 0];[1 1 0 0 0 0 0];[0 0 1 0 0 0 0];[0 0 1 1 0 0 0];[0 0 0 0 1 0 0];[0 0 0 0 0 1 1];[0 0 0 0 0 0 0]]

# ---- Initial state-space
x0 = [500, 4, 110, 300, 2, 20, 90]

# ---- Highest orders of reaction for tau-leaping
hor = [2,2,2,2,1,2,2]

# ---- Bounds on the initial search space for parameters
# ---- Real parameters are [.38, .04, .082, .12, .021, .1,.005, 13.21]
bounds = [[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], [1.0,1.0,1.0,1.0,1.0,1.0,1.0,20.0]]

# ---- Hazard functions
function F(x, θ, a, t)
	H(x, θ, a, t)
	a .*= θ
end

function H(x, θ, a, t)
	a[1] = 1.0
    a[2] = x[1]
    a[3] = x[1] * x[2]
    a[4] = x[3]
    a[5] = x[3] * x[4]
    a[6] = x[5]
    a[7] = x[6] * x[7]
    a[8] = 1.0
end
# ---- Configure system in SPICE
system = System(Model(x0,F,H,ν,bounds,hor=hor, obsname=[:R,:L,:RL,:G,:Ga,:Gbg,:Gd]), "./data/50-sets",routine=CEM(ssa=:Direct, nElite = 10, nRepeat = 1, nSamples=1000, maxIter=250, mSamples=20000, shoot=false, splitting=false, tauOpt=TauOpt(ϵ=0.1)))

# ---- Call parameter estimation routine
estimate(system, 1, "yeast")
