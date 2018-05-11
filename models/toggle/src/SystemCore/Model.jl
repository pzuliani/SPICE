struct Model
    x0::Vector{Int}
    F::Function
    H::Function
    ν::Matrix{Int}
    bounds::Vector{Vector{Float64}}
    mapping::Vector{Vector{Int}}
    obs::Vector{Int}
    obsname::Vector{Symbol}
    hor::Vector{Int}
    ns::Int
    nr::Int
    ps::Vector{Int}
    Model(x0::Vector{Int}, F::Function, H::Function, ν::Matrix{Int}, bounds::Vector{Vector{Float64}},mapping::Vector{Vector{Int}}, obs::Vector{Int}, obsname::Vector{Symbol}, hor::Vector{Int}, ns::Int, nr::Int, ps::Vector{Int}) = new(copy(x0), F, H, copy(ν), bounds, copy(mapping), copy(obs), obsname, copy(hor), ns, nr, ps)
end

Model(x0::Vector{Int}, F::Function, H::Function, ν::Matrix{Int}, bounds::Vector{Vector{Float64}}; obs::Vector{Int} = collect(1:length(x0)), obsname::Vector{Symbol} = Symbol[], hor::Vector{Int} = fill(2, length(x0)), ns::Int = length(x0), nr::Int = size(ν)[1], mapping::Vector{Vector{Int}} = [Int[i] for i in 1:length(bounds[1])], ps::Vector{Int} = collect(1:size(ν)[1])) = Model(convert(Vector{Int}, x0), F, H, convert(Matrix{Int}, ν), bounds, mapping, obs, obsname, hor, ns, nr, ps)

function Base.show(io::IO, m::Model)
    print("Model:\n\t\tx0: ", m.x0, ",\n\t\tObservables: ", m.obs, ",\n\t\tFunction: ", m.F, ",\n\t\tState Change Matrix: ", m.ν, ",\n\t\tHOR: ", m.hor, ",\n\t\t", m.ns, " species,\n\t\t", m.nr, " reactions.")
end