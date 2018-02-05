include("SSA.jl")

function addSim(sys::System, tauVar::TauVar, n::Int64)
    ens = newEnsemble(sys, n)  
    npoints = length(sys.data.dd)
    if sys.routine.ssa == :Direct
        for i in 2:npoints
            if i == 2 && sys.routine.shoot
                reShoot!(sys, ens, i-1)
            end
            _sim!(sys, ens, i)
            for k in eachindex(ens)
                for xx in ens[k].x
                    push!(ens[k].xa,xx)
                end
                push!(ens[k].ta,sys.times[i])
            end
        end
    elseif sys.routine.ssa == :Tau
        for i in 2:npoints
            if i == 2 && sys.routine.shoot
                reShoot!(sys, ens, i-1)
            end
            _sim!(sys, ens, i, tauVar)
            for k in eachindex(ens)
                for xx in ens[k].x
                    push!(ens[k].xa,xx)
                end
                push!(ens[k].ta,sys.times[i])
            end
        end
    end
    ens
end

function sim(sys::System, θ::Vector{Float64}, σ2::Vector{Float64}, n::Int64)
    sys.state.θ = θ
    sys.state.σ2 = diagm(σ2)
    tauVar = preAllocTauVar(sys)
    ens = Ensemble()
    append!(ens, addSim(sys, tauVar, n))
    ens
end
