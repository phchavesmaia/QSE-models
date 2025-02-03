function eq_plots(Lᵢ,wᵢ,rₙ,Pₙ;suffix::String="",dpi::Int64=900)
    " 
    This function is dependent on importing CairoMakie as Mke
    "
    
    ## plotting figure 2
    fig = Mke.Figure(size = (1100, 1000))
    # population
    ax = Mke.Axis(fig[1, 1], 
                title = "Log population",
                titlealign = :left, 
                xlabel = "Longitude", 
                ylabel = "Latitude")
    viz!(ax, fund.geometry, color = log.(Lᵢ), colormap = :viridis)
    Mke.Colorbar(fig[1, 2], limits = (minimum(log.(Lᵢ)), maximum(log.(Lᵢ))), 
            colormap = :viridis, flipaxis = true, ticks=range(round(minimum(log.(Lᵢ))), round(maximum(log.(Lᵢ))), length=5), flip_vertical_label=true)
    # wages
    ax = Mke.Axis(fig[1, 3], 
                title = "Log wages",
                titlealign = :left, 
                xlabel = "Longitude", 
                ylabel = "Latitude")
    viz!(ax, fund.geometry, color = log.(wᵢ), colormap = :viridis)
    Mke.Colorbar(fig[1, 4], limits = (minimum(log.(wᵢ)), maximum(log.(wᵢ))), 
            colormap = :viridis, flipaxis = true, ticks=range(round(minimum(log.(wᵢ))), round(maximum(log.(wᵢ))), length=5), label="Log points", flip_vertical_label=true)
    # land prices
    ax = Mke.Axis(fig[2, 1], 
                title = "Log land prices",
                titlealign = :left, 
                xlabel = "Longitude", 
                ylabel = "Latitude")
    viz!(ax, fund.geometry, color = log.(rₙ), colormap = :viridis)
    Mke.Colorbar(fig[2, 2], limits = (minimum(log.(rₙ)), maximum(log.(rₙ))), 
            colormap = :viridis, flipaxis = true, ticks=range(round(minimum(log.(rₙ))), round(maximum(log.(rₙ))), length=5), flip_vertical_label=true)
    # price index
    ax = Mke.Axis(fig[2, 3], 
                title = "Log price index",
                titlealign = :left, 
                xlabel = "Longitude", 
                ylabel = "Latitude")
    viz!(ax, fund.geometry, color = log.(Pₙ), colormap = :viridis)
    Mke.Colorbar(fig[2, 4], limits = (minimum(log.(Pₙ)), maximum(log.(Pₙ))), 
            colormap = :viridis, flipaxis = true, ticks=range(round(minimum(log.(Pₙ))), round(maximum(log.(Pₙ))), length=5), label="Log points", flip_vertical_label=true)
 
    return Mke.save(string("./figures/equilibrium",suffix,".png"), fig, px_per_unit = dpi/96)

end