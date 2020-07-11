function selectReaction(a::Vector{Float64}, suma::Float64, nr::Int)
    t = rand() * suma
    rxn = 0
    while t >= 0.0 && rxn < nr
        rxn += 1
        t -= a[rxn]
    end
    rxn
end

function selectReaction(a::SubArray, suma::Float64, nr::Int)
    t = rand() * suma
    rxn = 0
    while t >= 0.0 && rxn < nr
        rxn += 1
        t -= a[rxn]
    end
    rxn
end

function directMethod!(p::Path, F::Function, ν::Matrix{Int}, tf::Float64, nr::Int, mapping::Vector{Vector{Int}})
    #a = zeros(Float64, nr)
    fly = zeros(nr)
    while true
        # stop state-space exploration
        moleculeLimit!(p, tf) && break

        # Propensities
        F(p)
        suma = sum(p.a)

        if suma <= 0.0
            break
        end

        # Time step
        dt = randexp()/suma
        if (p.t + dt) > tf
            p.t = tf
            break
        end
        p.t += dt

        # calc on the fly
        @simd for i in eachindex(fly)
            @inbounds fly[i] += p.a[i]*dt
        end

        # Reaction selection & state update
        rxn = selectReaction(p.a, suma, nr)
        p.ra[rxn] += 1
        dx = view(ν, rxn, :)
        @simd for i in eachindex(p.x)
            @inbounds p.x[i] += dx[i]
        end

    end
    for j in eachindex(p.θ)
        for j2 in mapping[j]
            p.fly[j2] += fly[j2] ./ p.θ[j]
        end
    end
end
directMethod = directMethod!

function nDirectMethod!(p::Path, F::Function, ν::Matrix{Int}, tf::Float64, nr::Int, nd::Int, mapping::Vector{Vector{Int}})
    nstep = 0
    fly = zeros(nr)
    while nstep < nd
        moleculeLimit!(p, tf) && break

        # Propensities
        F(p)

        suma = sum(p.a)
        if suma <= 0.0
            break
        end

        # Time step
        dt = randexp()/suma
        if (p.t + dt) > tf
            p.t = tf
            break
        end
        p.t += dt

        # calc on the fly
        @simd for i in eachindex(fly)
            @inbounds fly[i] += p.a[i]*dt
        end
        # Reaction selection & state update
        rxn = selectReaction(p.a, suma, nr)
        p.ra[rxn] += 1
        dx = view(ν, rxn, :)
        @simd for i in eachindex(p.x)
            @inbounds p.x[i] += dx[i]
        end

        nstep += 1
    end
    for j in eachindex(p.θ)
        for j2 in mapping[j]
            p.fly[j2] += fly[j2] ./ p.θ[j]
        end
    end
end

function moleculeLimit!(p::Path, tf::Float64)
    if any(p.x .> 1500)
        p.t = tf
        p.dm .+= Inf
        return true
    end
    return false
end
