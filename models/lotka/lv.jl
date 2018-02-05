include("../../src/SPICE.jl")
using SPICE

# ---- State-change matrix
ν = [1 0; -1 1; 0 -1]

# ---- Initial state-space
x0 = [50,50]

# ---- Highest orders of reaction for tau-leaping
hor = [2,2]

# ---- Bounds on the initial search space for parameters
# ---- Real parameters are [0.5,0.0025,0.3]
bounds = [[0.0,0.0,0.0],[1.0,1.0,1.0]]

# ---- Hazard functions
function F(x, θ, a, t)
    H(x,θ,a,t)
    a .*= θ
end

function H(x, θ, a, t)
    a[1] = x[1]
    a[2] = x[1] * x[2]
    a[3] = x[2]
end


system = System(Model(x0,F,H, ν,bounds,hor=hor,obsname=[:X1,:X2]), "./data/5",routine=CEM(ssa=:Tau, nElite = 10, nRepeat=1, nSamples=1000, maxIter=250, mSamples=20000, shoot=true, splitting=false, tauOpt=TauOpt(ϵ=0.1)))
estimate(system, 1, "lotka")
