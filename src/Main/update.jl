# call this function to complete all the iteration update steps
function updateParameters!(sys::System, ens::Ensemble)
	# calculate useful vectors using CE-update formula
    rxnActual, rxnExpect = computeRxnStats(sys, ens)

    # update hyper-parameters
    meanHyper, varHyper = hyperUpdate(sys, ens)

    # update kinetic rate parameters
    θstep = computeRateMeans!(sys, rxnActual, rxnExpect, meanHyper)
    
    # update variances
    computeRateVar!(sys, ens, rxnActual, rxnExpect, θstep, varHyper)
    append!(sys.state.θs, sys.state.θ)
end

function computeRxnStats(sys::System, ens::Ensemble)
    len = length(ens)

    pl = length(sys.model.ps)
    # ---- pre-allocate
    rxnActual = zeros(Int64, (pl, len))
    rxnExpect = zeros(Float64, (pl, len))
    a  = zeros(Float64, sys.model.nr)

    for k in eachindex(ens)
        for myj in sys.model.ps
            rxnActual[myj, k] += sum(ens[k].ra[sys.model.mapping[myj]])
        end
        τvec = diff(ens[k].ta)
        for n in eachindex(τvec)
            sys.model.H(ens[k].xa[(n - 1) * sys.model.ns + 1:n*sys.model.ns], ens[k].θ, a, ens[k].ta[n])
            for myj in sys.model.ps
                rxnExpect[myj, k] += sum((τvec[n] .* a)[sys.model.mapping[myj]])
            end
        end
        sys.model.H(ens[k].x, ens[k].θ, a, ens[k].ta[end])
        for myj in sys.model.ps
            rxnExpect[myj, k] += sum((sys.times[end] - ens[k].t) .* a[sys.model.mapping[myj]])
        end
    end

    return rxnActual, rxnExpect
end

function hyperUpdate(sys::System, ens::Ensemble)
    n = length(sys.model.bounds[1]) - length(sys.model.ps)
    meanHyper = Vector{Float64}(n)
    varHyper = Vector{Float64}(n)

    for j in eachindex(meanHyper)
        v = Float64[]
        for k in eachindex(ens)
            push!(v,ens[k].θ[j+length(sys.model.ps)])
        end
        meanHyper[j] = mean(v)
        varHyper[j] = var(v)
    end
    return meanHyper, varHyper
end

function computeRateMeans!(sys::System, rxnActual::Matrix{Int}, rxnExpect::Matrix{Float64}, arm)
    θstep = vec(sum(rxnActual, 2) ./ sum(rxnExpect, 2))
    append!(θstep, arm)
    if sys.state.i == 1
        sys.state.θ = θstep
        return θstep
    else
        sys.state.θ .= sys.routine.smoothingParameter .* θstep .+ (1.0 - sys.routine.smoothingParameter) .* sys.state.θ
        return θstep
    end
end

function computeRateVar!(sys::System, ens::Ensemble, rxnActual::Matrix{Int}, rxnExpect::Matrix{Float64}, θstep::Vector{Float64}, ars)
    σ2step = vec(mean((rxnActual ./ max.(rxnExpect,1) .- θstep[sys.model.ps]).^2,2))
    append!(σ2step, ars)

    σ2step = max.(0.04.*θstep.^2, σ2step)
    β = 0.8 - 0.8 * (1.0 - 1.0 / sys.state.i)^5

    if sys.state.i == 1
        sys.state.σ2 = diagm(σ2step)
    else
        β = 0.8 - 0.8 * (1.0 - 1.0 / sys.state.i)^5
        sys.state.σ2 .= β .* diagm(σ2step) .+ (1.0 - β) .* sys.state.σ2
    end
end

function computeRateUncert!(sys::System, ens::Ensemble, rxnActual::Matrix{Int}, rxnExpect::Matrix{Float64}, θstep::Vector{Float64})
    len = length(ens)
    _term1 = diagm(zeros(sys.model.nr))
    _term2 = diagm(zeros(sys.model.nr))
    _pterm3 = zeros(sys.model.nr)
    for k in eachindex(ens)
        _term1 += diagm(rxnActual[:, k] ./ ens[k].θ.^2)
        _xterm = (rxnActual[:, k] ./ ens[k].θ) .- rxnExpect[:, k]
        _term2 += _xterm * _xterm'
        _pterm3 += _xterm
    end
    _term1 /= len
    _term2 /= len
    _pterm3 /= len
    _term3 = _pterm3 * _pterm3'
    s = _term1 + _term2 - _term3
    σ2 = inv(s)
    return σ2
end

