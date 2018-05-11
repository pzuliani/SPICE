include("dm.jl")
include("tau_aux.jl")
include("tau.jl")
function simIntDM!(sys::System, pathEnsemble::Ensemble, i::Int64)
    for k in eachindex(pathEnsemble)
        directMethod!(pathEnsemble[k], sys.model.F, sys.model.ν, sys.times[i], sys.model.nr)
    end
end

function simIntTau!(sys::System, pathEnsemble::Ensemble, i::Int64, tauVar::TauVar)
    for k in eachindex(pathEnsemble)
        optimisedTauLeap!(sys.model, pathEnsemble[k], sys.times[i], sys.routine.tauOpt, tauVar)
    end
end
