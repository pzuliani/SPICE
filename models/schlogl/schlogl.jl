include("../../src/SPICE.jl")
using .SPICE

# ---- State-change matrix
const ν = Array([1 -1 1 -1]')

# ---- Initial state-space
x0 = [250]

# ---- Highest orders of reaction for tau-leaping
const hor = [33]

# ---- Bounds on the initial search space for parameters
# ---- Real parameters are [3e-7,1e-4,1e-3,3.5]
const bounds = [[1e-9,1e-6,1e-5,1e-2],[1e-5, 1e-2, 1e-1, 100.0]]


function F(p)
	H(p)
	p.a .*= p.θ
end

# ---- Hazard functions
function H(p)
    p.a[1] = p.x[1] * (p.x[1]-1) / 2 * 100000
    p.a[2] = p.x[1] * (p.x[1]-1) * (p.x[1] - 2) / 6
    p.a[3] = 200000
    p.a[4] = p.x[1]
end


# how many times to run?
system = System(Model(x0,F,H, ν,bounds,hor=hor,obsname=[:X]), "./data/10",routine=CEM(ssa=:Tau, nElite = 10, nRepeat=10, nSamples=1000, maxIter=250, mSamples=20000, shoot=false, splitting=false, tauOpt=TauOpt(ϵ=0.1)))

estimate(system, 100, "schlogl")
