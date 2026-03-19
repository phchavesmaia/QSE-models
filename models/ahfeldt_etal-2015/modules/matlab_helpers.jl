module MatlabHelpers

export mapit, read_mat

using GeoStats, GeoIO, MAT
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

    # read data
    fileIn = matopen("$(path)/prepdata_big_TD$year.mat");
    dset = read(fileIn); close(fileIn);
    "
    It follows a brief data dictionary:
    ----------
    nobsXX = number of observations
    floorXX = rent prices
    empwplXX = workplace employment (population)
    emprsdXX = residential employment (population)
    ttXX = bilateral travel time matrix 
        τᵢⱼ is constructed on top of ttXX s.t. rows (i) denote residences and columns (j) denote workplaces.
    bzk86XX = mapping of Blocks to Bezirkes
    ----------
    "
    
    # defining suffix
    year=="86" ? suffix = "rw" : suffix="";
    
    # importing variables
    Qⱼ = vec(dset["floor$(year)$(suffix)"]); Hₘⱼ = vec(dset["empwpl$(year)$(suffix)"]) ; Hᵣᵢ = vec(dset["emprsd$(year)$(suffix)"]); # using the paper notation
    
    # distvar transposition (see https://github.com/Ahlfeldt/ARSW2015-toolkit/pull/2/changes/7068bddce63553eb9fb645d5f839ee4c17fbf9f9 for the argument for the transposition)  
    year=="86" ? τᵢⱼ = Matrix(dset["tt$(year)$(suffix)"]') : τᵢⱼ = dset["tt$(year)$(suffix)"];

    # area
    year=="86" ? Kᵢ = nothing : Kᵢ = vec(dset["area06"]);

    # block map
    block_bzk = vec(Int.(dset["bzk$(year)$(suffix)"])); # map from blocks to districts

    # close data
    dset = nothing;
    
    return Qⱼ, Hₘⱼ, Hᵣᵢ, τᵢⱼ, Kᵢ, block_bzk
end

end