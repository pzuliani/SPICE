include("./src/SPICE.jl")
using .SPICE

# ---- State-change matrix
const ν = [[1 0];[-1 0];[0 1];[0 -1]]
# ---- Initial state-space
x0 = [50,2500]

# ---- Highest orders of reaction for tau-leaping
const hor = [1,1]

# ---- Bounds on the initial search space for parameters
const bounds = [[1,1e-3,1,1e-3,50,50,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3], [50,1,50,1,1e4,500,500,10,10,10,10,10]]

# ---- Hazard functions
function F(p)
	H(p)
	p.a .*= p.θ[1:4]
end

function H(p)
    p.a[1] = p.θ[7] / (1 + p.θ[7] + p.θ[8]*p.x[2]^2)
    p.a[2] = p.x[1]
    p.a[3] = p.θ[10] / (1 + p.θ[10]+ p.θ[11]*p.x[1]^2)
    p.a[4] = p.x[2]
end

# ---- Configure system in SPICE
system = System(Model(x0,F,H,ν,bounds,hor=hor, ps=[1,2,3,4], obsname=[:X1, :X2]), "./data/iptg",routine=CEM(ssa=:Tau, nElite = 10, nRepeat = 1, nSamples=1000, maxIter=250, mSamples=20000, shoot=false, splitting=false, tauOpt=TauOpt(ϵ=0.1)))

# ---- Call parameter estimation routine
estimate(system, 1, "toggle")
