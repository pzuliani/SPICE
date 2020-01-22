function addFirstEns(sys::System, tauVar::TauVar)
    ens = addEnsemble(sys, method = :Sobol)
    npoints = length(sys.times)
    if sys.routine.ssa == :Direct
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

function addEns(sys::System, tauVar::TauVar, n::Int; gn::Int=0)
    ens = addEnsemble(sys, nSamples = n, method = :Dist, gn=gn)
    npoints = length(sys.times)
    if sys.routine.ssa == :Direct
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

function reShoot!(sys::System, ens::Ensemble, i::Int)
    for k in eachindex(ens)
        s = sample(1:length(sys.data), Weights(1.0./ens[k].dm))
        ens[k].x = vec(copy(sys.data[s][i,2:end]))
        ens[k].t = sys.times[i]
    end
end


function split!(sys::System, ens::Ensemble)
    inds = Vector{Int}(undef, 0)
    costs = Vector{Float64}(undef, length(ens))
    for i in eachindex(sys.data)
        for k in eachindex(ens)
            costs[k] = getCostVal(ens[k], i)
        end
        δ = quantile(costs, sys.routine.nElite / length(ens))
        append!(inds, findall(costs .<= δ))
    end
    s = sample(inds, length(ens))
    new_ens = Ensemble()
    ens = setPath(sys,ens,s)
end

function setPath(sys, ens, s)
    new_ens = Ensemble()
    if sys.state.i == 1
        sob = SobolSeq(length(sys.model.bounds[1]), sys.model.bounds[1],  sys.model.bounds[2])
        for i in eachindex(s)
            par = next!(sob)
            push!(new_ens, copyPath2(ens[s[i]], par))
        end
    else
        if sys.routine.sampling == :Normal
            dist = MvNormal(sys.state.θ, sys.state.σ2)
            for i in eachindex(s)
                par = rand(dist)
                push!(new_ens, copyPath2(ens[s[i]], par))
            end
        elseif sys.routine.sampling == :log
            dist = MvLogNormal(sys.state.θ, sys.state.σ2)
            for i in eachindex(s)
                par = rand(dist)
                push!(new_ens, copyPath2(ens[s[i]], par))
            end
        end
    end
    new_ens
end

function copyPath2(p::Path, par::Vector{Float64})
  	p2 = copyPath(p)
  	p2.θ = copy(par)
  	p2
end

function getEliteData(sys::System, ens::Ensemble)
    eliteSet = Ensemble()
    mincost = Inf
    δm = -Inf
    for i in eachindex(sys.data)
        costs = Vector{Float64}(undef, length(ens))
        for k in eachindex(ens)
            costs[k] = getCostVal(ens[k], i)
            costs[k] < mincost && (mincost = costs[k])
        end
        δ = quantile(costs, sys.routine.nElite / length(ens))
        if δ > δm
            δm = δ
        end
        for xx in ens[costs .<= δ]
            push!(eliteSet, xx)
        end
    end
    return eliteSet, δm, mincost
end

function getEliteGroups(sys::System, ens::Ensemble)
	ng = ens[end].group
	groupcosts = zeros(ng)
	for g in 1:ng
	    for i in eachindex(sys.data)
	    	mincost = Inf
	        for k in eachindex(ens)
	            (g == ens[k].group) && getCostVal(ens[k], i) < mincost && (mincost = getCostVal(ens[k], i))
	        end
	        groupcosts[g] += mincost
	    end
	end
	δ = quantile(groupcosts, sys.routine.nElite / ng)
	gelite = collect(1:ng)[groupcosts .<= δ]

	groupSet = Ensemble()

    for xx in ens
    	if in(xx.group, gelite)
        	push!(groupSet, xx)
        end
    end

    return groupSet, δ, minimum(groupcosts)
end
