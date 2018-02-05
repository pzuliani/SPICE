include("./src/SPICE.jl")
using SPICE

# ---- State-change matrix
ν = [[1 0];[-1 0];[0 1];[0 -1]]
# ---- Initial state-space
x0 = [50,2500]

# ---- Highest orders of reaction for tau-leaping
hor = [1,1]

# ---- Bounds on the initial search space for parameters
bounds = [[1.0,0,1.0,0,50,50,0,100], [100.0,1,100.0,1,500,500,100,1000]]

# ---- Hazard functions
function F(x, θ, a, t)
	H(Int64.(x), θ, a, t)
	a .*= θ[1:4]
end

function H(x, θ, a, t)
    a[1] = 1 / (1 + x[2]^2)
    a[2] = x[1]
    a[3] = 1 / (1 + x[1]^2)
    a[4] = (1+1) * x[2]
end

# ---- Configure system in SPICE
system = System(Model(x0,F,H,ν,bounds,hor=hor, ps=[1,2,3,4], obsname=[:X1, :X2]), "./data/iptg",routine=CEM(ssa=:Tau, nElite = 50, nRepeat = 1, nSamples=1000, maxIter=250, mSamples=20000, shoot=false, splitting=true, tauOpt=TauOpt(ϵ=0.15)))

# ---- Call parameter estimation routine
estimate(system, 1, "toggle-tau")
