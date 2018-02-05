immutable Model
    x0::Vector{Int32}
    F::Function
    H::Function
    ν::Matrix{Int32}
    bounds::Vector{Vector{Float64}}
    mapping::Vector{Vector{Int64}}
    obs::Vector{Int}
    obsname::Vector{Symbol}
    hor::Vector{Int64}
    ns::Int64
    nr::Int64
    ps::Vector{Int}
    Model(x0::Vector{Int32}, F::Function, H::Function, ν::Matrix{Int32}, bounds::Vector{Vector{Float64}},mapping::Vector{Vector{Int64}}, obs::Vector{Int64}, obsname::Vector{Symbol}, hor::Vector{Int64}, ns::Int64, nr::Int64, ps::Vector{Int}) = new(copy(x0), F, H, copy(ν), bounds, copy(mapping), copy(obs), obsname, copy(hor), ns, nr, ps)
end

Model(x0::Vector{Int64}, F::Function, H::Function, ν::Matrix{Int64}, bounds::Vector{Vector{Float64}}; obs::Vector{Int64} = collect(1:length(x0)), obsname::Vector{Symbol} = Symbol[], hor::Vector{Int64} = fill(2, length(x0)), ns::Int64 = length(x0), nr::Int64 = size(ν)[1], mapping::Vector{Vector{Int64}} = [Int64[i] for i in 1:length(bounds[1])], ps::Vector{Int} = collect(1:size(ν)[1])) = Model(convert(Vector{Int32}, x0), F, H, convert(Matrix{Int32}, ν), bounds, mapping, obs, obsname, hor, ns, nr, ps)

function Base.show(io::IO, m::Model)
    print("Model:\n\t\tx0: ", m.x0, ",\n\t\tObservables: ", m.obs, ",\n\t\tFunction: ", m.F, ",\n\t\tState Change Matrix: ", m.ν, ",\n\t\tHOR: ", m.hor, ",\n\t\t", m.ns, " species,\n\t\t", m.nr, " reactions.")
end