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

# ****************************
# **** Transaction Costs  **** 
# ****************************

# determining the graph representative of the 30x30 grid
N=30
diagonal_cost = 2^0.5;
straight_cost = 1.0;
A, A_weights = adjacency_matrix(N, diagonal_cost, straight_cost)
G = Graph(A)
gplot(G)

# estimating the least cost path matrix
for j ∈ 1:N*N
    dist = dijkstra_shortest_paths(G,j,A_weights).dists
    if j==1
        πᵢⱼ = dist
    else
        πᵢⱼ = hcat(πᵢⱼ,dist)
    end
end

# ***********************
# **** Productivity  **** 
# ***********************
using CSV ,DataFrames, GeoStats
import CairoMakie as Mke

Aᵢ = CSV.read("a.csv", DataFrame; header = false) 
rename!(Aᵢ,[:Productivity])

# plotting figure 1
grid = CartesianGrid(30, 30)
Aᵢ = georef(Aᵢ, grid)

fig = Mke.Figure()
viz(fig[1,1], Aᵢ.geometry, color = Aᵢ.Productivity, colormap = "viridis")
cbar(fig[1,2], Aᵢ.Productivity, colormap = "viridis")
fig
