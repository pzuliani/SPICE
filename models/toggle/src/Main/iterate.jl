function iterate!(sys::System, tauVar::TauVar)
    sys.state.i += 1
    ens = Ensemble()
    # ---- Is it the first iteration?
    if sys.state.i == 1
        cost = Inf
        mincost = Inf
        push!(sys.state.γs, cost)
        push!(sys.state.γmins, mincost)
        push!(sys.state.ns, 0)
        append!(ens, addFirstEns(sys, tauVar))

    # ---- All other iterations
    else
        append!(ens, addEns(sys, tauVar, round(Int,sys.state.ns[end]/sys.routine.nRepeat)))
    end
    
    nPaths = length(ens)

    while true
        nAdd = round(Int, min(nPaths/(3*sys.routine.nRepeat),(sys.routine.mSamples-nPaths)/sys.routine.nRepeat))

        # ---- Find & sort elite samples, & cost scores
        if sys.routine.nRepeat > 1
            eliteEnsemble, cost, mincost = getEliteGroups(sys, ens)
        else
            eliteEnsemble, cost, mincost = getEliteData(sys, ens)
        end

        # ---- Series of convergence checks
        if sys.state.i >= 3

            # ---- Has the cost function stagnated?
            if mincost == sys.state.γmins[end] == sys.state.γmins[end-1] 
                println("Terminating because S*(t)==S*(t-1)==S*(t-2)")
                
                updatestep!(sys, cost, mincost, eliteEnsemble, nPaths)
                return true

            # ---- Has the cost function stagnated?
            elseif cost == sys.state.γs[end] == sys.state.γs[end-1] 
                println("Terminating because γ(t)==γ(t-1)==γ(t-2)")

                updatestep!(sys, cost, mincost, eliteEnsemble, nPaths)
                return true

            # ---- Did we improve the cost function?
            elseif cost < minimum(sys.state.γs)
                updatestep!(sys, cost, mincost, eliteEnsemble, nPaths)
                return false

            # ---- Did we improve the best sample?
            elseif mincost < minimum(sys.state.γmins)
                updatestep!(sys, cost, mincost, eliteEnsemble, nPaths)
                return false

            # ---- Are we consistently hitting the sample limit?
            elseif sys.routine.mSamples <= nPaths == sys.state.ns[end] == sys.state.ns[end-1] 
                println("Terminating because N(t)==N(t-1)==N(t-2)")

                nullstep!(sys, cost, mincost, eliteEnsemble, nPaths)
                return true

            # ---- Do we need more samples?
            elseif nPaths < sys.routine.mSamples
                gn = ens[end].group
                nPaths += nAdd * sys.routine.nRepeat
                append!(ens, addEns(sys, tauVar, nAdd, gn=gn))
                continue
        
            # ---- Did we hit the max sample limit?
            elseif nPaths >= sys.routine.mSamples
                nullstep!(sys, cost, mincost, eliteEnsemble, nPaths)
                return false
            end
        else

            # ---- Did we improve the cost function?
            if cost < minimum(sys.state.γs)
                #update
                updatestep!(sys, cost, mincost, eliteEnsemble, nPaths)
                return false

            # ---- Did we improve the best sample?
            elseif mincost < minimum(sys.state.γmins)
                updatestep!(sys, cost, mincost, eliteEnsemble, nPaths)
                return false

            # ---- Do we need more samples?
            elseif nPaths < sys.routine.mSamples
                gn = ens[end].group
                nPaths += nAdd * sys.routine.nRepeat
                append!(ens, addEns(sys, tauVar, nAdd, gn=gn))
                continue
            
            # ---- Did we hit the max sample limit?
            elseif nPaths >= sys.routine.mSamples
                nullstep!(sys, cost, mincost, eliteEnsemble, nPaths)
                return false
            end
        end
    end
end

function updatestep!(sys::System, cost::Float64, mincost::Float64, eliteEnsemble::Ensemble, nIt::Int)
    push!(sys.state.γs, cost)
    push!(sys.state.γmins, mincost)
    updateParameters!(sys, eliteEnsemble)
    push!(sys.state.ns, nIt)
    println(sys.state)
end
function nullstep!(sys::System, cost::Float64, mincost::Float64, eliteEnsemble::Ensemble, nIt::Int)
    push!(sys.state.γs, cost)
    push!(sys.state.γmins, mincost)
    append!(sys.state.θs, sys.state.θ)
    push!(sys.state.ns, nIt)
    println(sys.state)
end