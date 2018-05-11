function _optimisedTauLeap!(model::Model, p::Path, tf::Float32, tauOpt::TauOpt, tauVar::TauVar)
    
    # ---- Main loop
    a = zeros(Float64, model.nr)
    while p.t < tf
        moleculeLimit!(p, tf) && break

        model.F(p.x, p.θ, a, p.t)
        if sum(a) <= 0.0
            break
        end
        tr = tf - p.t

        τ, switch = computeLeap(model, p, a, tr, tauOpt, tauVar)
        
        # ---- Switch to direct method?
        if switch == true
            _nDirectMethod!(p, model.F, model.ν, tf, model.nr, tauOpt.nd)

        # ---- Else accept step
        else
            for i in eachindex(p.x)
                p.x[i] += tauVar.dx[i]
            end
            p.t += τ
        end
    end
end

_optimisedTauLeap = _optimisedTauLeap!
