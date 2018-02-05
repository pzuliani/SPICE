function iterate!(sys::System, tauVar::TauVar)
    sys.state.i += 1
    ens = Ensemble()
    nAdd = round(Int,sys.routine.nSamples/3)
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

        # ---- Find & sort elite samples, & cost scores
        eliteEnsemble, cost, mincost = getEliteData(sys, ens)
        ens = nothing

        # ---- Series of convergence checks
        if sys.state.i >= 3

            # ---- Has the cost function stagnated?
            if mincost == sys.state.γmins[end] == sys.state.γmins[end-1] # == sys.state.γmins[end-2] == sys.state.γmins[end-3]
                println("Terminating because S*(t)==S*(t-1)==S*(t-2)")
                
                updatestep!(sys, cost, mincost, eliteEnsemble, nPaths)
                return true

            # ---- Has the cost function stagnated?
            elseif cost == sys.state.γs[end] == sys.state.γs[end-1] # == sys.state.γs[end-2] == sys.state.γs[end-3]
                println("Terminating because γ(t)==γ(t-1)==γ(t-2)")

                updatestep!(sys, cost, mincost, eliteEnsemble, nPaths)
                return true

            # ---- Are we consistently hitting the sample limit?
            elseif sys.routine.mSamples <= nPaths == sys.state.ns[end] == sys.state.ns[end-1] # == sys.state.ns[end-2] == sys.state.ns[end-3]
                println("Terminating because N(t)==N(t-1)==N(t-2)")

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

            # ---- Do we need more samples?
            elseif nPaths < sys.routine.mSamples
                ens = eliteEnsemble
                nPaths += nAdd * sys.routine.nRepeat
                append!(ens, addEns(sys, tauVar, nAdd))
                continue
        
            # ---- Did we hit the max sample limit?
            elseif nPaths >= sys.routine.mSamples
                updatestep!(sys, cost, mincost, eliteEnsemble, nPaths)
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
                ens = eliteEnsemble
                nPaths += nAdd * sys.routine.nRepeat
                append!(ens, addEns(sys, tauVar, nAdd))
                continue
            
            # ---- Did we hit the max sample limit?
            elseif nPaths >= sys.routine.mSamples
                updatestep!(sys, cost, mincost, eliteEnsemble, nPaths)
                return false
            end
        end
    end
end

function updatestep!(sys::System, cost::Float64, mincost::Float64, eliteEnsemble::Ensemble, nIt::Int64)
    push!(sys.state.γs, cost)
    push!(sys.state.γmins, mincost)
    updateParameters!(sys, eliteEnsemble)
    push!(sys.state.ns, nIt)
    println(sys.state)
end