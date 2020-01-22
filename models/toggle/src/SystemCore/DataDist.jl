struct DataDistribution
    dd::Dict{Float64, Distributions.Distribution}
end

# assuming normally distributed data
function createDataDistribution(dfs::Vector{DataFrame}, obs::Vector{Symbol})
    nobs = length(obs)

    # extract unique timepoints
    times = Float64[]
    for df in dfs
        append!(times, convert(Vector{Float64}, df[:Time]))
    end
    times = sort!(unique(times))

    # create vector of data distributions over time
    npoints = length(times)
    dd = Dict{Float64, Distributions.Distribution}()
    for t in times
        ar = Float64[]
        for df in dfs
            if any(df[:Time] .== t)
                ar = vcat(ar, convert(Array{Float64}, df[df[:Time] .== t, obs]))
            end
        end
        m = vec(mean(ar, 1))
        s = vec(std(ar, 1))

        # if std = NaN? Implement 10% std as a rough approximation
        if any(s .== NaN)
            s[s .== NaN] = m[s .== NaN] ./ 10
        end
        if any(s .== 0.0)
            s[s .== 0.0] = m[s .== 0.0] ./ 10
        end

        dd[Float64(t)] = MvNormal(m, s)
    end
    return DataDistribution(dd), convert(Vector{Float64},times)
end

function createDataDistributionM(dfs::Vector{DataFrame}, obs::Vector{Symbol})
    nobs = length(obs)

    # extract unique timepoints
    times = Float64[]
    for df in dfs
        append!(times, convert(Vector{Float64}, df[!, :Time]))
    end
    times2 = sort!(unique(times))

    # create vector of data distributions over time
    npoints = length(times2)
    dd = Dict{Float64, Distributions.Distribution}()
    for t in times2
        ar = MvNormal[]
        for df in dfs
            if any(df[!, :Time] .== t)
                for aa in 1:size(df[df[!, :Time] .== t, obs])[1]
                    vv = vec(convert(Array{Float64},df[df[!, :Time] .== t, obs][aa,:]))
                    push!(ar,MvNormal(vv, vv./10))
                end
            end
        end
        dd[Float64(t)] = MixtureModel(vec(ar))
    end
    return DataDistribution(dd), times2
end

function genstate(d::DataDistribution, t::Float64)
    xa = round.(Int,rand(d.dd[t], 1))
    while any(xa.<0)
        xa = round.(Int,rand(d.dd[t], 1))
    end
    vec(xa)
end

function loadData(path::String)
    datasets = DataFrame[] # Init empty vector
    if ispath(path)
        if isfile(path)
            push!(datasets, CSV.read(path))
            return datasets
        else
            files = readdir(path) # grab all files in the given folder
            for f in files
                push!(datasets, CSV.read(path * "/" * f))
            end
            return datasets
        end
    else
        error("Enter a valid file or folder")
    end
end
