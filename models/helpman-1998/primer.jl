using Graphs, GraphPlot, LinearAlgebra, SparseArrays

# *********************
# **** Load Files  **** 
# *********************

function load_dir(dir::String)
    files = readdir(dir)
    for file in files
        if endswith(file, ".jl")
            include(joinpath(dir, file))
        end
    end
end

cd("D:/Dropbox/learn-julia/qse/models/helpman-1998/")
load_dir("functions")

# ********************
# **** Distances  **** 
# ********************

# determining the graph representative of the 30x30 grid
N=30
diagonal_cost = 2^0.5
straight_cost = 1.0
A, A_weights = adjacency_matrix(N, diagonal_cost, straight_cost)
G = Graph(A)
gplot(G)

# estimating the "least cost path (distance)" matrix 
distₙᵢ = zeros(N*N,N*N)
for j ∈ 1:N*N
    distₙᵢ[:,j] = dijkstra_shortest_paths(G,j,A_weights).dists
end

# ***********************
# **** Productivity  **** 
# ***********************
using CSV ,DataFrames, GeoStats, CairoMakie

Aᵢ = CSV.read("a.csv", DataFrame; header = false) 
rename!(Aᵢ,[:Productivity])
grid = CartesianGrid(30, 30)
Aᵢ = georef(Aᵢ, grid)

# plotting figure 1
fig = Figure()
ax = Axis(fig[1, 1], 
              title = "Log productivity",
              titlealign = :left, 
              xlabel = "Longitude", 
              ylabel = "Latitude")
viz!(ax, Aᵢ.geometry, color = Aᵢ.Productivity, colormap = :viridis)
Colorbar(fig[1, 2], limits = (minimum(Aᵢ.Productivity), maximum(Aᵢ.Productivity)), 
        colormap = :viridis, flipaxis = true, ticks=-3:1:3, label="Log points", flip_vertical_label=true)
fig

# *********************
# **** Parameters  **** 
# *********************

α = 0.25 # Davis & Ortalo-Magne (2011)
σ = 5 # Simonovska & Waugh (2014)
ϕ = 0.375 # Based on international trade data
dₙᵢ = distₙᵢ .^ ϕ # trade costs are a constant elasticity function of effective distance by assumption 
τⁱⁿ = 2 # by assumption

Aᵢ = @transform(Aᵢ, :center = centroid(:geometry))
@transform(Aᵢ, :X = first(coordinates(:center)))