@compat abstract type Routine end #cross-version compatability

struct TauOpt
    nc::Int64
    ϵ::Float64
    dtf::Int64
    nd::Int64

    TauOpt(; nc::Int64 = 10, ϵ::Float64 = 0.1, dtf::Int64 = 10, nd::Int64 = 100) = new(nc, ϵ, dtf, nd)
end

struct CEM <: Routine
    #Routine parameters
    nSamples::Int64
    nRepeat::Int64
    nElite::Int64
    smoothingParameter::Float64
    shoot::Bool
    splitting::Bool
    sampling::Symbol

    #SSA + Parameters
    ssa::Symbol
    tauOpt::TauOpt

    #Termination parameters
    mSamples::Int64
    maxIter::Int64

    function CEM(; nSamples::Int64 = 1000, nRepeat::Int64 = 1, nElite::Int64 = 10, smoothingParameter::Float64 = 0.7, shoot = false, splitting = false, sampling = :log, ssa::Symbol = :Direct, tauOpt::TauOpt = TauOpt(), mSamples::Int64 = nSamples*50, maxIter::Int64 = 250)
        new(nSamples, nRepeat, nElite, smoothingParameter, shoot, splitting, sampling, ssa, tauOpt, mSamples, maxIter)
    end
end

function Base.show(io::IO, routine::CEM)
    print("\n\t\tnSamples: ", routine.nSamples, ",\n\t\tsmoothing: ", routine.smoothingParameter, ",\n\t\tSSA: ", routine.ssa, ".")
end