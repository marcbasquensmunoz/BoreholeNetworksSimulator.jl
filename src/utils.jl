
function get_all_julia_files_in_dir(dir)
    files_or_dirs = readdir(dir; join=true)
    return filter(Base.contains(r".jl$"), vcat([isdir(f) ? get_all_julia_files_in_dir(f) : f for f in files_or_dirs]...))
end

function sort_dependencies!(d, order)
    function sorting(s) 
        index = findall(i -> Base.contains(s, i), order)
        if length(index) >= 1
            return index[1]
        else
            return length(d)+1
        end
    end
    sort!(d, by=sorting)
end