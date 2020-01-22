const Ensemble = Vector{Path}

function Base.show(io::IO, ens::Ensemble)
    println(length(ens), " Trajectories")
end

function addEnsemble(sys::System; nSamples::Int = sys.routine.nSamples, method::Symbol = :Sobol, gn = 0)
    ens = Ensemble(undef, nSamples*sys.routine.nRepeat)
    if method == :Sobol
        if sys.routine.sampling == :Normal
            if sys.state.i == 1
                sob = SobolSeq(length(sys.model.bounds[1]), sys.model.bounds[1],  sys.model.bounds[2])
            else
                low = sys.state.θ .- 2 .* sqrt.(diag(sys.state.σ2))
                up = sys.state.θ .+ 2 .* sqrt.(diag(sys.state.σ2))
                sob = SobolSeq(length(sys.model.bounds[1]), low, up)
            end
            setEnsemble!(sys, ens, sob, gn)
        elseif sys.routine.sampling == :log
            if sys.state.i == 1
                sob = SobolSeq(length(sys.model.bounds[1]), log.(sys.model.bounds[1]),  log.(sys.model.bounds[2]))
            else
                low = sys.state.θ .- 2 .* sqrt.(diag(sys.state.σ2))
                up = sys.state.θ .+ 2 .* sqrt.(diag(sys.state.σ2))
                sob = SobolSeq(length(sys.model.bounds[1]), low, up)
            end
            setEnsembleLog!(sys, ens, sob, gn)
        end
    elseif method == :Dist
        if sys.routine.sampling == :Normal
            dist = MvNormal(sys.state.θ, sys.state.σ2)
            setEnsemble!(sys, ens, dist, gn)
        elseif sys.routine.sampling == :log
            dist = MvLogNormal(sys.state.θ, sys.state.σ2)
            setEnsemble!(sys, ens, dist, gn)
        end
    end
    return ens
end

function setEnsemble!(sys::System, ens::Ensemble, sob::Sobol.ScaledSobolSeq, gn)
    l = round(Int,length(ens)/sys.routine.nRepeat,RoundUp)
    j = 0
    for i in gn+1:gn+l
        ps = next!(sob)
        for r in 1:sys.routine.nRepeat
            j+=1
            if j > length(ens)
                return
            end
            ens[j] = Path(i, ps, [0,0], sys.times[1], zeros(Int, sys.model.nr))
        end
    end
end

function setEnsemble!(sys::System, ens::Ensemble, dist, gn)
    l = round(Int,length(ens)/sys.routine.nRepeat,RoundUp)
    j = 0
    for i in gn+1:gn+l
        ps = rand(dist)
        for r in 1:sys.routine.nRepeat
            sample = genstate(sys.data, sys.times[1])
            for s in eachindex(sample)
                sample[s] = round(Int32,sample[s]/ps[4+s])
            end
            j+=1
            if j > length(ens)
                return
            end
            ens[j] = Path(i, ps, sample, sys.times[1], zeros(Int, sys.model.nr))
        end
    end
end


function setEnsembleLog!(sys::System, ens::Ensemble, sob::Sobol.ScaledSobolSeq, gn)
    l = round(Int,length(ens)/sys.routine.nRepeat,RoundUp)
    j = 0
    for i in gn+1:gn+l
        ps = exp.(next!(sob))
        for r in 1:sys.routine.nRepeat
            sample = genstate(sys.data, sys.times[1])
            for s in eachindex(sample)
                sample[s] = round(Int32,sample[s]/ps[4+s])
            end
            j+=1
            if j > length(ens)
                return
            end
            ens[j] = Path(i, ps, sample, sys.times[1], zeros(Int, sys.model.nr))
        end
    end
end
