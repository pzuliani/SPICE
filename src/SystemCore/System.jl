include("Model.jl")
include("Routine.jl")
include("SystemState.jl")
include("DataDist.jl")


mutable struct System
    model::Model
    routine::CEM
    data::Vector{Array}
    times::Vector{Float64}
    state::SystemState

    function System(model::Model, folder::String; routine::CEM = CEM(), state::SystemState = SystemState()) 
        ds = loadData(folder)
        da = arData(ds)
        times = Float64.(da[1][:,1])
        new(model, routine, da, times, state)
    end
end

function Base.show(io::IO, s::System)
    print("\nSystem:", "\n\t", s.model, "\n\tRoutine: ", s.routine, "\n\tState: ", s.state)
end

function resetSystem!(sys::System)
    sys.state = SystemState()
end

include("Path.jl")
include("Ensemble.jl")
