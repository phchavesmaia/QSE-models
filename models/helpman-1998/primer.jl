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

try 
    cd("D:/Dropbox/learn-julia/qse/models/helpman-1998/")
catch
    cd("C:/Users/pedro.maia/Dropbox/learn-julia/qse/models/helpman-1998/")
end
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
using CSV, DataFrames, GeoStats, CairoMakie

# read productivity data
data = CSV.read("a.csv", DataFrame; header = false)
rename!(data,[:Aᵢ])
# indicate cells identifiers
data.cell = 1:size(data,1) 
# construct cells geometries
grid = CartesianGrid(30, 30)
data = georef(data, grid)
# indicate country flags
aux = values(data)
aux.Country = data.cell.∈ Ref(geosplit(data,0.5,(1,0))[1].cell)
data = georef(aux, grid)

# plotting figure 1
fig = Figure()
ax = Axis(fig[1, 1], 
              title = "Log productivity",
              titlealign = :left, 
              xlabel = "Longitude", 
              ylabel = "Latitude")
viz!(ax, data.geometry, color = data.Aᵢ, colormap = :viridis)
Colorbar(fig[1, 2], limits = (minimum(data.Aᵢ), maximum(data.Aᵢ)), 
        colormap = :viridis, flipaxis = true, ticks=-3:1:3, label="Log points", flip_vertical_label=true)
fig

# *********************
# **** Parameters  **** 
# *********************

α = 0.25 # Davis & Ortalo-Magne (2011)
σ = 5.0 # Simonovska & Waugh (2014)
ϕ = 0.375 # Based on international trade data
dₙᵢ = distₙᵢ .^ ϕ # trade costs are a constant elasticity function of effective distance by assumption 
τᵒ = τⁱ = 2 # by assumption


