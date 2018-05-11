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
bounds = [[1e-6,1e-6,1e-6],[10.0,10.0,10.0]]

# ---- Hazard functions
function F(p)
    H(p)
    p.a .*= p.θ
end

function H(p)
    p.a[1] = p.x[1]
    p.a[2] = p.x[1] * p.x[2]
    p.a[3] = p.x[2]
end


system = System(Model(x0,F,H, ν,bounds,hor=hor,obsname=[:X1],obs=[1]), "./data/5",routine=CEM(ssa=:Direct, nElite = 10, nRepeat=1, nSamples=1000, maxIter=250, mSamples=20000, shoot=false, splitting=false, sampling=:log, tauOpt=TauOpt(ϵ=0.1)))
estimate(system, 100, "lotka")
