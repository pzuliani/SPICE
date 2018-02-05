@compat abstract type Routine end #cross-version compatability

immutable TauOpt
    nc::Int64
    ϵ::Float64
    dtf::Int64
    nd::Int64

    TauOpt(; nc::Int64 = 10, ϵ::Float64 = 0.2, dtf::Int64 = 10, nd::Int64 = 100) = new(nc, ϵ, dtf, nd)
end

immutable CEM <: Routine
    #Routine parameters
    nSamples::Int64
    nRepeat::Int64
    nElite::Int64
    smoothingParameter::Float64
    shoot::Bool
    splitting::Bool

    #SSA + Parameters
    ssa::Symbol
    tauOpt::TauOpt

    #Termination parameters
    mSamples::Int64
    maxIter::Int64

    function CEM(; nSamples::Int64 = 500, nRepeat::Int64 = 10, nElite::Int64 = 50, smoothingParameter::Float64 = 0.9, shoot = true, splitting = false, ssa::Symbol = :Direct, tauOpt::TauOpt = TauOpt(), mSamples::Int64 = nSamples*100, maxIter::Int64 = 500)
        new(nSamples, nRepeat, nElite, smoothingParameter, shoot, splitting, ssa, tauOpt, mSamples, maxIter)
    end
end

function Base.show(io::IO, routine::CEM)
    print("\n\t\tnSamples: ", routine.nSamples, ",\n\t\tsmoothing: ", routine.smoothingParameter, ",\n\t\tSSA: ", routine.ssa, ".")
end