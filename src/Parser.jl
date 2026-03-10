import Base.+

"""
    +(buffer::Array{String,1}, line::String)

Extend `Base.+` to append a `String` to a `String` array in-place.
This allows using `+(buffer, line)` as a concise alternative to `push!(buffer, line)`.

### Arguments
- `buffer::Array{String,1}`: the string array to append to.
- `line::String`: the string to append.
"""
function +(buffer::Array{String,1}, line::String)
    push!(buffer, line)
end

"""
    read_reaction_file(path_to_file::String) -> Array{String,1}

Read in aVFF reaction file and return an array of records as strings.
The comments are lines starting with "//".

### Arguments
- `path_to_file::String`: the path to the VFF reaction file.

### Returns
- `Array{String,1}`: an array of records as strings.
"""
function read_reaction_file(path_to_file::String)::Array{String,1}

    # initialize -
    vff_file_buffer = String[]
    vff_reaction_array = Array{String,1}()

    # Read in the file -
    open("$(path_to_file)", "r") do file
        for line in eachline(file)
            +(vff_file_buffer, line)
        end
    end

    # process -
    for reaction_line ∈ vff_file_buffer
        
        # skip comments and empty lines -
        if (occursin("//", reaction_line) == false && 
            isempty(reaction_line) == false)
        
            # grab -
            push!(vff_reaction_array,reaction_line)
        end
    end

    # return -
    return vff_reaction_array
end