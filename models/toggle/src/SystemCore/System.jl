include("Model.jl")
include("Routine.jl")
include("SystemState.jl")
include("DataDist.jl")


type System
    model::Model
    routine::CEM
    data::DataDistribution
    times::Vector{Float64}
    state::SystemState

    function System(model::Model, folder::String; routine::CEM = CEM(), state::SystemState = SystemState()) 
        create = createDataDistributionM(loadData(folder), model.obsname)
        new(model, routine, create[1], create[2], state)
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
