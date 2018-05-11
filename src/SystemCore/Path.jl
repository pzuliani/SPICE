mutable struct Path
	group::Int
    θ::Vector{Float64}
    x::Vector{Int}
    t::Float64
    a::Vector{Float64}
    fly::Vector{Float64}
    ra::Vector{Int}
    dm::Vector{Float64}
    w::Float64
    Path(group::Int,θ::Vector{Float64}, x::Vector{Int}, t::Float64, a::Vector{Float64}, fly::Vector{Float64}, ra::Vector{Int}, dm::Vector{Float64}, w::Float64) = new(group, θ, x, t, a, fly, ra, dm, w)
end

Path(group::Int,θ::Vector{Float64}, x::Vector{Int}, t::Float64, ra::Vector{Int}, nd::Int) = Path(group, θ, copy(x), t, zeros(Float64, length(ra)), Float64.(ra), ra, zeros(nd), 1.0)
Path(group::Int,θ::Vector{Float64}, x::Vector{Int}, t::Float64, ra::Vector{Int}, dm::Vector{Float64}) = Path(group, θ, copy(x), t, zeros(Float64, length(ra)), Float64.(ra), ra, dm, 1.0)

function Base.show(io::IO, p::Path)
    print("Path { \n group: ", p.group, " \n θ: ", p.θ, " \n x: ", p.x, ", \n t: ", p.t,  ", \n reaction firings: ", p.ra, ", \n objective value: ", p.dm, "}")
end

function getCostVal(p::Path)
    minimum(p.dm)
end
function getCostVal(p::Path, i)
    p.dm[i]
end

function getRates(p::Path)
    p.θ
end

function getGroup(p::Path)
	p.group
end

function getGroupCounts(ens::Vector{Path})
	counts(getGroup.(ens))
end

function overwritePath!(p1::Path, p2::Path)
	p1.group = p2.group
    p1.θ = copy(p2.θ)
    p1.x = copy(p2.x)
    p1.t = p2.t
    p1.a = copy(p2.a)
    p1.fly = copy(p2.fly)
    p1.ra = copy(p2.ra)
    p1.dm = copy(p2.dm)
end

function copyPath(p::Path)
	group = p.group
    θ = copy(p.θ)
    x = copy(p.x)
    t = p.t
    a = copy(p.a)
    fly = copy(p.fly)
    ra = copy(p.ra)
    dm = copy(p.dm)
    w = copy(p.w)
    Path(group, θ, x, t, a, fly, ra, dm, w)
end
function copyPath(p::Path, rw::Int)
	group = p.group
    θ = copy(p.θ)
    x = copy(p.x)
    t = p.t
    a = copy(p.a)
    fly = copy(p.fly)
    ra = copy(p.ra)
    dm = copy(p.dm)
    w = copy(p.w) / rw
    Path(group, θ, x, t, a, fly, ra, dm, w)
end