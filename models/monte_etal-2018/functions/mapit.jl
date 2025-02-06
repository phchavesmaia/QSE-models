function mapit(path_from,data,title; path_to::String="", dpi::Int64=900, label_legend::String="Log points")
    
    # import geometries
    shapes = GeoIO.load(path_from)

    # make figure
    fig = Mke.Figure(size=(500,600))
    ax = Mke.Axis(fig[1, 1], 
                title = title,
                titlealign = :left)
    viz!(ax, shapes.geometry, color = vec(data), colormap = :viridis, showsegments=true, segmentsize=.2 )
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