# call this function to complete all the iteration update steps
function updateParameters!(sys::System, ens::Ensemble)
	# calculate useful vectors using CE-update formula
    rxnActual, rxnExpect = computeRxnStats(sys, ens)

    # update kinetic rate parameters
    θstep = computeRateMeans!(sys, rxnActual, rxnExpect)

    # update variances
    computeRateVar!(sys, ens, rxnActual, rxnExpect, θstep)
    append!(sys.state.θs, sys.state.θ)
end

function iszeroreaction(p::Path)
    any(p.ra .== 0)
end

function computeRxnStats(sys::System, ens::Ensemble)
    len = length(ens)

    pl = length(sys.model.ps)
    # ---- pre-allocate
    rxnActual = zeros(Float64, (pl, len))
    rxnExpect = zeros(Float64, (pl, len))
    a  = zeros(Float64, sys.model.nr)

    for k in eachindex(ens)
        sys.model.F(ens[k])
        for j in sys.model.ps
            @inbounds for j2 in sys.model.mapping[j]
                rxnActual[j, k] += ens[k].ra[j2]
                rxnExpect[j, k] += ens[k].fly[j2]
                rxnExpect[j, k] += (sys.times[end] - ens[k].t) .*  ens[k].a[j2]
            end
        end

    end

    return rxnActual, rxnExpect
end

function computeRateMeans!(sys::System, rxnActual::Matrix{Float64}, rxnExpect::Matrix{Float64})
    θstep = vec(sum(rxnActual, dims = 2) ./ sum(rxnExpect, dims = 2))

    if sys.routine.sampling == :log
        θstep = log.(θstep)
    end
    if sys.state.i == 1
        sys.state.θ = θstep
        return θstep
    else
        sys.state.θ .= sys.routine.smoothingParameter .* θstep .+ (1.0 - sys.routine.smoothingParameter) .* sys.state.θ
        return θstep
    end
end

function computeRateVar!(sys::System, ens::Ensemble, rxnActual::Matrix{Float64}, rxnExpect::Matrix{Float64}, θstep::Vector{Float64})
    if sys.routine.sampling == :log
        tmp = log.((rxnActual .+ 0.01) ./ rxnExpect)
        σ2step = vec(var(tmp, dims = 2))
    else
        σ2step = vec(var(tmp, dims = 2))
    end

    if sys.state.i == 1
        sys.state.σ2 = diagm(σ2step)
    else
        β = 0.8 - 0.8 * (1.0 - 1.0 / sys.state.i)^5
        sys.state.σ2 .= β .* diagm(σ2step) .+ (1.0 - β) .* sys.state.σ2
    end
end
