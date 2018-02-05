include("../../src/SPICE.jl")
using SPICE

# ---- State-change matrix
ν = [1 -1 1 -1]'

# ---- Initial state-space
x0 = [250]

# ---- Highest orders of reaction for tau-leaping
hor = [33]

# ---- Bounds on the initial search space for parameters
# ---- Real parameters are [3e-7,1e-4,1e-3,3.5]
bounds = [[1e-8,1e-5,1e-4,1e-1],[1e-6, 1e-3, 1e-2, 10.0]]


function F(x, θ, a, t)
	H(x, a)
	a .*= θ
end

# ---- Hazard functions
function H(x, θ, a, t)
    a[1] = x[1] * (x[1]-1) / 2 * 100000 
    a[2] = x[1] * (x[1]-1) * (x[1] - 2) / 6
    a[3] = 200000
    a[4] = x[1]
end


# how many times to run? 
system = System(Model(x0,F,H, ν,bounds,hor=hor,obsname=[:X]), "./data/10",routine=CEM(ssa=:Tau, nElite = 10, nRepeat=1, nSamples=1000, maxIter=250, mSamples=20000, shoot=true, splitting=false, tauOpt=TauOpt(ϵ=0.1)))

estimate(system, 50, "schlogl-gen")
