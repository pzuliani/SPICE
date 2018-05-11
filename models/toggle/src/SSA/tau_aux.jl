type TauVar
    neg_nu::Matrix{Int}
    g::Vector{Float64}
    irs::Vector{Int}
    dx::Vector{Int}
    l::Vector{Float64}
    k::Vector{Int}

    TauVar(neg_nu::Matrix{Int}, g::Vector{Float64}, irs::Vector{Int}, dx::Vector{Int}, l::Vector{Int}, k::Vector{Int}) = new(neg_nu, g, irs, dx, l, k)
end

function preAllocTauVar(sys::System)
    # useful precomputations
    neg_nu = copy(sys.model.ν)
    neg_nu[sys.model.ν .>= zero(eltype(sys.model.ν))] = zero(eltype(neg_nu))
    
    g = zeros(Float64, sys.model.ns)
    for i in eachindex(sys.model.hor)
        sys.model.hor[i] == 0 && (g[i] = 0.0; continue)
        sys.model.hor[i] == 1 && (g[i] = 1.0; continue)
        sys.model.hor[i] == 2 && (g[i] = 2.0; continue)
    end

    irsBits = vec(any(sys.model.ν .!= zero(eltype(sys.model.ν)), 1))
    irs = Vector{Int}(0)
    for i in eachindex(irsBits)
        irsBits[i] && push!(irs, i)
    end

    dx = zeros(eltype(sys.model.x0), sys.model.ns)
    l = zeros(Int, sys.model.nr)
    k = zeros(Int, sys.model.nr)

    TauVar(neg_nu, g, irs, dx, l, k)
end