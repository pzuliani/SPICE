include("../../src/SPICE.jl")
using .SPICE

# SIR Model
# 0 -> S    (m*(S+I+R))
# S -> 0    (m*S)
# I -> 0      ((m+v)*I)
# R -> 0    (m*R)
# S+I -> 2I (b*S*I)
# I -> R     (r*I)

# ---- State-change matrix
const ν = [1 0 0; # 0 -> S
        -1 0 0; # S -> 0
        0 -1 0; # I -> 0 (m * I)
        0 -1 0; # I -> 0 (v * I)
        0 0 -1; # R -> 0
        -1 1 0; # S+I -> 2I
        0 -1 1] # I -> R

# ---- Initial state-space
# S, I, R
x0 = [1000, 10, 0]

# ---- Highest orders of reaction for tau-leaping
# hor(S) = 2 due to reaction propensity [b*S*I]
# hor(I) = 2 due to reaction propensity [b*S*I]
# hor(R) = 1 due to reaction propensity [m*R]
const hor = [2, 2, 1]

# ---- Bounds on the initial search space for parameters
const bounds = [[1e-6,1e-6,1e-6,1e-6],  # set to your lower bounds
                [1e-3,1e-3,1e-3,1e-3]]  # set to your upper bounds

# ---- Hazard/propensity functions
function F(p)
    H(p)
    p.a[1] *= p.θ[1] # m * ...
    p.a[2] *= p.θ[1] # m * ...
    p.a[3] *= p.θ[1] # m * ...
    p.a[4] *= p.θ[2] # v * ...
    p.a[5] *= p.θ[1] # m * ...
    p.a[6] *= p.θ[3] # b * ...
    p.a[7] *= p.θ[4] # r * ...
end

function H(p)
    p.a[1] = p.x[1] + p.x[2] + p.x[3] # ... (S + I + R)
    p.a[2] = p.x[1] # ... S
    p.a[3] = p.x[2] # ... I
    p.a[4] = p.x[2] # ... I
    p.a[5] = p.x[3] # ... R
    p.a[6] = p.x[1] * p.x[2] # ... S * I
    p.a[7] = p.x[2] # ... I
end

mapping = [[1, 2, 3, 5], [4], [6], [7]]

model = Model(x0, F, H, ν, bounds, hor = hor, obsname = [:S, :I, :R], mapping = mapping)

system = System(model, "./data", routine=CEM(ssa=:Direct, nElite = 10, nRepeat=1, nSamples=1000, maxIter=250, mSamples=20000, shoot=false, splitting=false, sampling=:log, tauOpt=TauOpt(ϵ=0.1)))

estimate(system, 1, "SIR")
