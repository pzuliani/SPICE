module SPICE

    # ---- Modules required
    using Distributions
    using DataFrames
    using Sobol
    using Compat
    using StatsBase
    using Statistics
    using LinearAlgebra
    using Random
    using CSV

    # ---- Imports & includes
    import Base.show, Base.copy
    include("SystemCore/System.jl")
    include("SSA/SSA.jl")
    include("Main/Main.jl")
    include("Sim/sim.jl")

    # ---- Export types
    export CEM,
    Model,
    DataSets,
    TauOpt,
    System,
    SystemState,

    # ---- Export functions
    est,
    estimate,
    resetSystem!,
    sim
end
