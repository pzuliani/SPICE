function est(sys::System)
    println(sys)
    println("Running...")
    term = false
    tauVar = preAllocTauVar(sys)

    # ---- Run while we haven't terminated!
    while !term && sys.state.i < sys.routine.maxIter
        term = iterate!(sys, tauVar)
    end
    println("Finished")
    return sys.state
end

function estimate(system::System, nrun::Int, output::String)
    for i in 1:nrun
        # ---- Where to save?
        filepathEst = "./results/"*output*"-estimates.csv"
        filepathTrace = "./results/traces/"*output*"$i-trace.csv"
        mkpath("./results/traces/")

        # ---- Start timer
        tic()
        state = est(system)

        # ---- End timer
        t = toc()

        # ---- Estimates
        f1 = open(filepathEst, "a")
        writecsv(f1, [state.θ' diag(state.σ2)' t])
        close(f1)

        # ---- Costs
        γs = state.γs[2:end]
        l = length(γs)

        # ---- Estimate traces
        θs = reshape(state.θs, (system.model.nr, l))'

        # ---- nSamples
        ns = state.ns[2:end]

        it = collect(1:l)
        _run = repmat([i],l)
        traceAr = [θs γs ns it _run]

        f2 = open(filepathTrace, "w")
        writecsv(f2, traceAr)
        close(f2)

        resetSystem!(system)
    end
    return true
end