function loadData(path::String)
    datasets = DataFrame[] # Init empty vector
    if ispath(path)
        if isfile(path)
            push!(datasets, readtable(path))
            return datasets
        else
            files = readdir(path) # grab all files in the given folder
            for f in files
                push!(datasets, readtable(path * "/" * f))
            end
            return datasets
        end
    else
        error("Enter a valid file or folder")
    end
end

function arData(datasets::Vector{DataFrame})
    return Array.(datasets)
end

