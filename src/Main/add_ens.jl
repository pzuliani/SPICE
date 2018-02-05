function addFirstEns(sys::System, tauVar::TauVar)
    ens = newFirstEnsemble(sys) 
    npoints = length(sys.data.dd)
    if sys.routine.ssa == :Direct
        #getCost!(sys, ens, 1)
        for i in 2:npoints
            if sys.routine.shoot #|| i == 2
                reShoot!(sys, ens, i-1)
            end
            if i > 2 && sys.routine.splitting
                split!(sys, ens)
            end
            simIntDM!(sys, ens, i)
            getCost!(sys, ens, i)
        end
    elseif sys.routine.ssa == :Tau
        #getCost!(sys, ens, 1)
        for i in 2:npoints
            if sys.routine.shoot #|| i == 2
                reShoot!(sys, ens, i-1)
            end
            if i > 2 && sys.routine.splitting
                split!(sys, ens)
            end
            simIntTau!(sys, ens, i, tauVar)
            getCost!(sys, ens, i)
        end
    end
    ens
end

function addEns(sys::System, tauVar::TauVar, n::Int64)
    ens = newEnsemble(sys, n)  
    npoints = length(sys.data.dd)
    if sys.routine.ssa == :Direct
        #getCost!(sys, ens, 1)
        for i in 2:npoints
            if sys.routine.shoot #|| i == 2
                reShoot!(sys, ens, i-1)
            end
            if i > 2 && sys.routine.splitting
                split!(sys, ens)
            end
            simIntDM!(sys, ens, i)
            getCost!(sys, ens, i)
        end
    elseif sys.routine.ssa == :Tau
        #getCost!(sys, ens, 1)
        for i in 2:npoints
            if sys.routine.shoot #|| i == 2
                reShoot!(sys, ens, i-1)
            end
            if i > 2 && sys.routine.splitting
                split!(sys, ens)
            end
            simIntTau!(sys, ens, i, tauVar)
            getCost!(sys, ens, i)
        end
    end
    ens
end

function reShoot!(sys::System, ens::Ensemble, i::Int64)
    for k in eachindex(ens)
        sample = genstate(sys.data, sys.times[i])
        for ii in eachindex(sample) 
            ens[k].x[ii] = sample[ii]
        end
        for xx in ens[k].x
            push!(ens[k].xa, xx)
        end
        ens[k].t = sys.times[i]
        push!(ens[k].ta, ens[k].t)
    end
end

function split!(sys::System, ens::Ensemble)
    inds = Vector{Int64}(0)
    costs = Vector{Float64}(length(ens))
    for k in eachindex(ens)
        costs[k] = getCostVal(ens[k])
    end
    δ = quantile(costs, (sys.routine.nElite / (sys.routine.nSamples*sys.routine.nRepeat)))
    append!(inds, find(costs .<= δ))
    cc = counts(inds, length(ens))
    for k in eachindex(ens)
        if cc[k] == 0
            i = sample(inds)
            overwritePath!(ens[k], ens[i])
        end
    end
end

function getEliteData(sys::System, ens::Ensemble)
    eliteSet = Ensemble()
    costs = Vector{Float64}(length(ens))
    mincost = Inf
    for k in eachindex(ens)
        costs[k] = getCostVal(ens[k])
        costs[k] < mincost && (mincost = costs[k])
    end
    δ = quantile(costs, (sys.routine.nElite / length(ens)))
    for xx in ens[costs .<= δ]
        push!(eliteSet, xx)
    end
    return eliteSet, δ, mincost
end
