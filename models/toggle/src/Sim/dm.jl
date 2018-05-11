function _directMethod!(path::Path, F::Function, ν::Matrix{Int32}, tf::Float32, nr::Int64)
    a = zeros(Float64, nr)
    while path.t < tf
        # stop silly state-space exploration
        moleculeLimit!(path, tf) && break

        # Propensities
        F(path.x, path.θ, a, path.t)
        suma = sum(a)
        if suma <= 0.0
            break
        end

        # Time step
        dt = Float32(rand(Exponential(1.0 / suma)))
        if path.t + dt > tf
            break
        end
        path.t += dt

        # Reaction selection & state update
        rxn = selectReaction(a, suma, nr)
        path.ra[rxn] += 1
        dx = view(ν, rxn, :)
        @simd for i in eachindex(path.x)
            @inbounds path.x[i] += dx[i]
        end

    end
end
_directMethod = _directMethod!

function _nDirectMethod!(p::Path, F::Function, ν::Matrix{Int32}, tf::Float32, nr::Int64, nd::Int64)
    nstep = 0
    a = zeros(Float64, nr)
    while p.t < tf && nstep < nd
        # Propensities

        F(p.x, p.θ, a, p.t)
        suma = sum(a)
        if suma <= 0.0
            break
        end

        # Time step
        dt = Float32(rand(Exponential(1.0 / suma)))
        p.t += dt

        # Reaction selection & state update
        rxn = selectReaction(a, suma, nr)
        p.ra[rxn] += 1
        dx = view(ν, rxn, :)
        @simd for i in eachindex(p.x)
            @inbounds p.x[i] += dx[i]
        end

        nstep += 1
    end
end