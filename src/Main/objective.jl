function getCost!(sys::System, pathEnsemble::Ensemble, i::Int)
    for k in eachindex(pathEnsemble)
        cost!(sys, pathEnsemble[k], i)
    end
end

function cost!(sys::System, p::Path, i::Int)
    for ds in eachindex(sys.data)
        for m in sys.model.obs
            @inbounds p.dm[ds] += ((sys.data[ds][i,m+1]-p.x[m]))^2
        end
    end
end