const Ensemble = Vector{Path}

function Base.show(io::IO, ens::Ensemble)
    println(length(ens), " Trajectories")
end

function newFirstEnsemble(sys::System)
    sob = SobolSeq(length(sys.model.bounds[1]), sys.model.bounds[1],  sys.model.bounds[2])
    ens = Ensemble(sys.routine.nSamples*sys.routine.nRepeat)
    setFirstEnsemble!(sys, ens, sob)
    return ens
end

function newFirstEnsembleLog(sys::System)
    sob = SobolSeq(length(sys.model.bounds[1]), log10.(sys.model.bounds[1]),  log10.(sys.model.bounds[2]))
    ens = Ensemble(sys.routine.nSamples*sys.routine.nRepeat)
    setFirstEnsembleLog!(sys, ens, sob)
    return ens
end

function setFirstEnsemble!(sys::System, ens::Ensemble, sob::Sobol.ScaledSobolSeq)
    for i in eachindex(ens)
        ps = next(sob)
        for r in 1:sys.routine.nRepeat
            ens[i] = Path(ps, sys.model.x0, sys.times[1], zeros(Int64, sys.model.nr), length(sys.data.dd))
        end
    end
end

function setFirstEnsembleLog!(sys::System, ens::Ensemble, sob::Sobol.ScaledSobolSeq)
    for i in eachindex(ens)
        ps = 10.^next(sob)
        for r in 1:sys.routine.nRepeat
            ens[i] = Path(ps, sys.model.x0, sys.times[1], zeros(Int64, sys.model.nr), length(sys.data.dd))
        end
    end
end

function newEnsemble(sys::System, n::Int64)
    dist = MvNormal(sys.state.θ, sys.state.σ2)
    ens = Ensemble(n*sys.routine.nRepeat)
    setEnsemble!(sys, ens, dist)
    return ens
end

function sampleDist(dist)
    v = rand(dist)
    while any(v .< 0)
        v = rand(dist)
    end
    return v
end

function setEnsemble!(sys::System, ens::Ensemble, dist)
    for i in eachindex(ens)
        ps = sampleDist(dist)
        for r in 1:sys.routine.nRepeat
            ens[i] = Path(ps, sys.model.x0, sys.times[1], zeros(Int64, sys.model.nr), length(sys.data.dd))
        end
    end
end
