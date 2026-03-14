module MatlabHelpers

export mapit, read_mat

using GeoStats, GeoIO
import CairoMakie as Mke

function mapit(path_from::String, data::Vector{Float64}, title::String; path_to::String="", dpi::Int64=900, label_legend::String="Log points")
    
    # import geometries
    shapes = GeoIO.load(path_from)

    # make figure
    fig = Mke.Figure()
    ax = Mke.Axis(fig[1, 1], 
                aspect=Mke.DataAspect(),
                title = title,
                titlealign = :left)
    viz!(ax, shapes.geometry, color = data, colormap = :viridis, showsegments=true, segmentsize=.2 )
    Mke.Colorbar(fig[1, 2], limits = (minimum(data), maximum(data)), 
        colormap = :viridis, flipaxis = true, label=label_legend, flip_vertical_label=true)
    Mke.hidedecorations!(ax)  # hides ticks, grid and lables
    
    if path_to == ""
        return fig
    else 
        Mke.save(path_to, fig, px_per_unit = dpi/96) #900 dpi
        return fig
    end
end

function read_mat(year::String, path::String="./data/input")
    fileIn = matopen("$(path)/prepdata_big_TD$year.mat");
    dset = read(fileIn); close(fileIn);
    "
    It follows a brief data dictionary:
    ----------
    nobsXX = number of observations
    floorXX = rent prices
    empwplXX = workplace employment (population)
    emprsdXXrw = residential employment (population)
    tt86rXX = bilateral travel time matrix s.t. rows (i) denote workplaces and columns (j) denote residences; which I transpose.
    bzk86XX = mapping of Blocks to Bezirkes
    ----------
    "
    Qⱼ = vec(dset["floor$(year)rw"]); Hₘⱼ = vec(dset["empwpl$(year)rw"]) ; Hᵣᵢ = vec(dset["emprsd$(year)rw"]); τᵢⱼ = Matrix(dset["tt$(year)rw"]'); # using the paper notation
    block_bzk = vec(Int.(dset["bzk$(year)rw"])); # map from blocks to districts 
    dset = nothing;
    return Qⱼ, Hₘⱼ, Hᵣᵢ, τᵢⱼ, block_bzk
end

end