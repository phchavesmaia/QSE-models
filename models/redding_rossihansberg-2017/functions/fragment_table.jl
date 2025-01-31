function fragment_table(labs, tab, file_name; path="./")
    "
    Recieves an matrix and returns the respective fragment table in LaTeX.
    Requires DelimitedFiles package.
    "
    if size(tab, 1) != size(labs, 1)
        throw("Row lables must be equal in size as the columns entries")
    elseif typeof(labs) != Vector{String}
        throw("Row lables must be delivered as a column vector")
    else
        tab_wrapped = ["&\$" * string(tab[i, j]) * "\$" for i in axes(tab,1), j in axes(tab,2)]
        ending = ["\\\\\\addlinespace" for i in axes(tab,1)]
        f = [labs tab_wrapped ending]
        open(string(path, file_name, ".tex"), "w") do io
            writedlm(io, f, " ")
        end
    end
end