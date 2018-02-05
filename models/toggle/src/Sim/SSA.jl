include("dm.jl")
include("tau.jl")

function _sim!(sys::System, pathEnsemble::Ensemble, i::Int64)
    for k in eachindex(pathEnsemble)
        _directMethod!(pathEnsemble[k], sys.model.F, sys.model.ν, sys.times[i], sys.model.nr)
    end
end

function _sim!(sys::System, pathEnsemble::Ensemble, i::Int64, tauVar::TauVar)
    for k in eachindex(pathEnsemble)
        _optimisedTauLeap!(sys.model, pathEnsemble[k], sys.times[i], sys.routine.tauOpt, tauVar)
    end
end
