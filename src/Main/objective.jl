function loglikelihood(dists, p::Path, consts::Vector{Float64})
    like = 0.0
    xx = reshape(p.x, (length(p.x), 1))
    for i in eachindex(dists)
        like += pdf(dists[i], xx)[1] * consts[i]
    end
    log(like)
end

function getCost!(sys::System, pathEnsemble::Ensemble, i::Int64)
    dists = components(sys.data.dd[sys.times[i]])
    n = ncomponents(sys.data.dd[sys.times[i]])
    consts = (2*pi).^(length.(dists)/2).*sqrt.(det.(diagm.(var.(dists))))./n
    for k in eachindex(pathEnsemble)
        costll!(dists, pathEnsemble[k], consts)
    end
end

function costll!(dists, p::Path, consts::Vector{Float64})
    p.dm -= loglikelihood(dists, p, consts)
end