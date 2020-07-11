function optimisedTauLeap!(model::Model, p::Path, tf::Float64, tauOpt::TauOpt, tauVar::TauVar)

    # ---- Main loop
    fly = zeros(model.nr)
    while p.t < tf
        moleculeLimit!(p, tf) && break
        model.F(p)
        if sum(p.a) <= 0.0
            break
        end
        tr = tf - p.t

        τ, switch = computeLeap(model, p, tr, tauOpt, tauVar)

        # ---- Switch to direct method?
        if switch == true
            nDirectMethod!(p, model.F, model.ν, tf, model.nr, tauOpt.nd, model.mapping)

        # ---- Else accept step
        else
            for i in eachindex(p.x)
                @inbounds p.x[i] += tauVar.dx[i]
            end

            @simd for i in eachindex(fly)
                @inbounds fly[i] += p.a[i]*τ
            end
            p.t += τ
        end
    end
    for j in eachindex(p.θ)
        for j2 in model.mapping[j]
            p.fly[j2] += fly[j2] ./ p.θ[j]
        end
    end
end

optimisedTauLeap = optimisedTauLeap!

function gUpdate!(hor::Vector{Int}, g::Vector{Float64}, x::Vector{Int})
    for i in eachindex(hor)
        if hor[i] == 22
            g[i] = 2.0 + 1.0 / (x[i] - 1.0)
            continue
        elseif hor[i] == 32
            g[i] = 1.5 * (2.0 + 1.0 / (x[i] - 1.0))
            continue
        elseif hor[i] == 33
            g[i] = 3.0 + 1.0 / (x[i] - 1.0) + 2.0 / (x[i] - 2.0)
        end
    end
end

function computeLeap(model::Model, p::Path, tr::Float64, tauOpt::TauOpt, tauVar::TauVar)::Tuple{Float64,Bool}
    τ = -Inf

    # ---- index of non-critical reactions
    Jncr = Vector{Int}(undef, 0)
    Jcr = Vector{Int}(undef, 0)

    for j in eachindex(tauVar.l)
        tauVar.l[j] = minimum(p.x ./ abs.(vec(tauVar.neg_nu[j, :])))
        tauVar.l[j] >= tauOpt.nc ? push!(Jncr, j) : push!(Jcr, j)
    end

    Acr = view(p.a, Jcr)
    Ancr = view(p.a, Jncr)

    if length(Jncr) > 0

        # ---- No need to update hor[i] == 1, hor[i] == 2 as is constant
        gUpdate!(model.hor, tauVar.g, p.x)
        μ = zeros(Float64, length(p.x))
        σ2 = zeros(Float64, length(p.x))
        for i in tauVar.irs
            @inbounds for j in Jncr
                μ[i] += model.ν[j, i] * p.a[j]
                σ2[i] += model.ν[j, i]^2 * p.a[j]
            end
        end

        maxTerm = max(maximum(view(tauOpt.ϵ .* (p.x ./ tauVar.g), tauVar.irs)), 1.0)
        leftTerm = (maxTerm ./ abs.(μ))
        rightTerm = (maxTerm.^2 ./ σ2)
        τ1 = minimum([leftTerm; rightTerm])
    else
        τ1 = Inf
    end

    switchCondition = (tauOpt.dtf / sum(p.a))::Float64

    # ---- total propensity of critical reactions
    sumAcr = sum(Acr)::Float64

    # ---- Loop until suitable τ that leads to non-negative states
    while true

        # ---- Check if we should switch to direct method
        if τ1 < switchCondition
            return τ, true
        end

        # ---- Time to next critical reaction
        τ2 = randexp()/sumAcr#rand(Exponential(1.0 / sumAcr))#-log(rand()) / sumAcr

        # ---- Case 1
        if τ1 <= τ2
            τ = τ1
            τ > tr && (τ = tr)
            @inbounds for s in eachindex(tauVar.k)
                tauVar.k[s] = zero(eltype(tauVar.k))
            end
            poissonSampler!(Jncr, tauVar.k, Ancr, τ)

        # ---- Case 2
        else
            τ = τ2
            τ > tr && (τ = tr)
            rxn = selectReaction(Acr, sumAcr, length(Acr))::Int
            @inbounds for s in eachindex(tauVar.k)
                tauVar.k[s] = zero(eltype(tauVar.k))
            end
            view(tauVar.k, Jcr)[rxn] = 1
            poissonSampler!(Jncr, tauVar.k, Ancr, τ)
        end

        # ---- Check proposed step, recompute τ if negative species, or break
        for i in eachindex(tauVar.dx)
            tauVar.dx[i] = 0
            for j in eachindex(tauVar.k)
                tauVar.dx[i] += tauVar.k[j] * model.ν[j, i]
            end
        end
        jump = true
        for i in eachindex(p.x)
            if p.x[i] + tauVar.dx[i] < zero(eltype(p.x))
                τ1 /= 2.0
                jump = false
                break
            end
        end
        if jump
            @simd for i in eachindex(p.ra)
                @inbounds p.ra[i] += tauVar.k[i]
            end
            break
        end
    end
    return τ, false
end

function poissonSampler!(Jncr::Vector{Int}, k::Vector{Int}, Ancr::SubArray, τ::Float64)
    for i in eachindex(Ancr)
        view(k, Jncr)[i] = Int(rand(Poisson(Ancr[i] * τ)))
    end
end
