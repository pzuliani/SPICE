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
bounds = [[1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6], [10.0,10.0,10.0,10.0,10.0,10.0,10.0,100.0]]

# ---- Hazard functions
function F(p)
	H(p)
	p.a .*= p.θ
end

function H(p)
	p.a[1] = 1.0
    p.a[2] = p.x[1]
    p.a[3] = p.x[1] * p.x[2]
    p.a[4] = p.x[3]
    p.a[5] = p.x[3] * p.x[4]
    p.a[6] = p.x[5]
    p.a[7] = p.x[6] * p.x[7]
    p.a[8] = 1.0
end
# ---- Configure system in SPICE
system = System(Model(x0,F,H,ν,bounds,hor=hor, obsname=[:R,:L,:RL,:G,:Ga,:Gbg,:Gd]), "./data/5-sets",routine=CEM(ssa=:Direct, nElite = 10, nRepeat = 1, nSamples=1000, maxIter=250, mSamples=20000, shoot=false, splitting=false, sampling=:log))

estimate(system, 100, "yeast")


