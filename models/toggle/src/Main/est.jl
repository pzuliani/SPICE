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

function estimate(system::System, nrun::Int64, output::String)
    for i in 1:nrun
        # ---- Where to save?
        filepathEst = "./results/"*output*"-estimates.csv"
        filepathTrace = "./results/traces/"*output*"$i-trace.csv"
        mkpath("./results/traces/")

        # ---- Start timer
        start_time = time()
        state = est(system)

        # ---- End timer
        stop_time = time()
        t = stop_time - start_time

        # ---- Estimates
        f1 = open(filepathEst, "a")
        CSV.write(f1, DataFrame([Array(state.θ') Array(diag(state.σ2)') t]), delim=',')
        close(f1)

        # ---- Costs
        γs = state.γs[2:end]
        l = length(γs)

        # ---- Estimate traces
        θs = Array(reshape(state.θs, (length(state.θ), l))')

        # ---- nSamples
        ns = state.ns[2:end]

        it = collect(1:l)
        _run = repeat([i],l)
        traceAr = [θs γs ns it _run]

        f2 = open(filepathTrace, "w")
        CSV.write(f2, DataFrame(traceAr), delim=',')
        close(f2)

        resetSystem!(system)
    end
    return true
end
