type SystemState
    i::Int64
    θ::Vector{Float64}
    σ2::Matrix{Float64}
    γs::Vector{Float64}
    γmins::Vector{Float64}
    ns::Vector{Int64}
    θs::Vector{Float64}
    
    SystemState(; i::Int64 = 0, θ::Vector{Float64} = Vector{Float64}(), σ2::Matrix{Float64} = Matrix{Float64}(0, 0), γs::Vector{Float64} = Vector{Float64}(), δmin_history::Vector{Float64} = Vector{Float64}(), ns::Vector{Float64} = Vector{Float64}(), θs::Vector{Float64} = Vector{Float64}()) = new(i, θ, σ2, γs, δmin_history, ns, θs)
end

function Base.show(io::IO, s::SystemState)
    if s.i == 0
        print("\n\t\tIterations: ", s.i, ",\n\t\tCurrent θ: ", s.θ, ",\n\t\tVariance ", diag(s.σ2), ",\n\t\tCost Value: Inf,\n\t\tMin Cost Value: Inf.")
    else
        print("\n\t\tIterations: ", s.i, ",\n\t\tCurrent θ: ", s.θ, ",\n\t\tVariance ", diag(s.σ2), ",\n\t\tCost Value: ", s.γs[end], ",\n\t\tMin Cost Value: ", s.γmins[end], ",\n\t\tIteration samples: ", s.ns[end], ",\n\t\tCumulative samples: ", sum(s.ns), ".")
    end
end