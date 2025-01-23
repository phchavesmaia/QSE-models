using Graphs, GraphPlot, LinearAlgebra, SparseArrays, StatsBase

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
dist = zeros(N*N,N*N)
for j ∈ 1:N*N
    dist[:,j] = dijkstra_shortest_paths(G,j,A_weights).dists
end

# ***********************
# **** Fundamentals  **** 
# ***********************
using CSV, DataFrames, GeoStats, CairoMakie

# read productivity data
fund = CSV.read("a.csv", DataFrame; header = false)
rename!(fund,[:Aᵢ])
fund.Aᵢ = exp.(fund.Aᵢ)
# indicate cells identifiers
fund.cell = 1:size(fund,1) 
# construct cells geometries
grid = CartesianGrid(30, 30)
fund = georef(fund, grid)
# indicate country flags
aux = values(fund)
aux.Country = fund.cell.∈ Ref(geosplit(fund,0.5,(1,0))[1].cell)
# cell area 
aux.Area = 100 * ones(N*N)
# back with the geometry
fund = georef(aux, grid)


# plotting figure 1
fig = Figure()
ax = Axis(fig[1, 1], 
              title = "Log productivity",
              titlealign = :left, 
              xlabel = "Longitude", 
              ylabel = "Latitude")
viz!(ax, fund.geometry, color = log.(fund.Aᵢ), colormap = :viridis)
Colorbar(fig[1, 2], limits = (minimum(log.(fund.Aᵢ)), maximum(log.(fund.Aᵢ))), 
        colormap = :viridis, flipaxis = true, ticks=-3:1:3, label="Log points", flip_vertical_label=true)
fig

# *********************
# **** Parameters  **** 
# *********************

α = 0.75 # Davis & Ortalo-Magne (2011)
σ = 5.0 # Simonovska & Waugh (2014)
ϕ = 0.375 # Based on international trade data
d = dist .^ ϕ # trade costs are a constant elasticity function of effective distance by assumption
d[diagind(d)] .= 1 # iceberg transport costs are one
# cell trade tax
τⁱ = ones(N*N,N*N).*2 # by assumption
τⁱ[diagind(τⁱ)] .= 1 # by construction
τⁱ = τⁱ .- 1
# country trade tax
τᵒ = ones(N*N,N*N)
τᵒ[fund.Country.==1,fund.Country.==0] .=2 # east to west
τᵒ[fund.Country.==0,fund.Country.==1] .=2 # west to east
τᵒ = τᵒ .- 1
F = 1.0 # production fixed cost
LL=153889 # US civilian labor force 2010 (Statistical Abstract, millions)
LLwest = LL*(sum(fund.Country.==1)/size(values(fund),1))
LLeast = LL*(sum(fund.Country.==0)/size(values(fund),1))

# ****************************************************************************
# **** Open Economy Solve for Endogenous Variables in Initial Equilibrium ****
# ****************************************************************************

println(">>>> Start Wage and Population Convergence <<<<")

wᵢ, Lᵢ, πₙᵢ, dom_πₙᵢ, converge, xtic = solveHLwCtyOpen_E(fund,d,τⁱ,τᵒ,N*N)

