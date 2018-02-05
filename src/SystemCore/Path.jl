type Path
    θ::Vector{Float64}
    x::Vector{Int32}
    xa::Vector{Int32}
    t::Float32
    ta::Vector{Float32}
    ra::Vector{Int64}
    dm::Float64

    Path(θ::Vector{Float64}, x::Vector{Int32}, xa::Vector{Int32}, t::Float32, ta::Vector{Float32}, ra::Vector{Int64}, dm::Float64) = new(θ, x, xa, t, ta, ra, dm)
end

Path(θ::Vector{Float64}, x::Vector{Int32}, t::Float32, ra::Vector{Int64}, i::Int) = Path(θ, copy(x), copy(x), Float32(t), Float32[t], ra, 0.0)

Base.copy(p::Path) = Path(p.θ, copy(p.x), copy(p.xa), p.t, copy(p.ta), copy(p.ra), copy(p.dm))

function Base.show(io::IO, p::Path)
    print("Path { \n θ: ", p.θ, " \n x: ", p.x, ", \n xa: ", p.xa, ", \n t: ", p.t, ", \n ta: ", p.ta, ", \n reaction firings: ", p.ra, ", \n objective value: ", p.dm, "}")
end

function getCostVal(p::Path)
    p.dm
end

function getRates(p::Path)
    p.θ
end

function overwritePath!(p1::Path, p2::Path)
    p1.x = copy(p2.x)
    p1.xa = copy(p2.xa)
    p1.t = copy(p2.t)
    p1.ta = copy(p2.ta)
    p1.ra = copy(p2.ra)
    p1.dm = copy(p2.dm)
end

# function splitPath!(p1::Path, p2::Path)
#     dist = MvNormal(p2.θ, 0.2*p2.θ)
#     p1.θ = sampleDist(dist)
#     p1.x = copy(p2.x)
#     p1.xa = copy(p2.xa)
#     p1.t = copy(p2.t)
#     p1.ta = copy(p2.ta)
#     p1.ra = copy(p2.ra)
#     p1.dm = copy(p2.dm)
# end

