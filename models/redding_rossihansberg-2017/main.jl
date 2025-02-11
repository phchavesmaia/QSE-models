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
    cd("D:/Dropbox/learn-julia/qse/models/redding_rossihansberg-2017/")
catch
    cd("C:/Users/pedro.maia/Dropbox/learn-julia/qse/models/redding_rossihansberg-2017/")
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
using CSV, DataFrames, GeoStats
import CairoMakie as Mke

# read productivity data 
fund = CSV.read("./data/a.csv", DataFrame; header = false);
rename!(fund,[:Aᵢ]);
fund.Aᵢ = vec(reshape(transpose(reshape(fund.Aᵢ,N,N)),N*N,1)); # matlab reshapes things in a different direction, fixing it so that the harmonization is correct
# indicate cells identifiers
fund.cell = 1:size(fund,1) ;
# construct cells geometries
grid = CartesianGrid(30, 30);
fund = georef(fund, grid);
# indicate country flags
aux = values(fund);
aux.Country = fund.cell.∈ Ref(geosplit(fund,0.5,(1,0))[1].cell);
# cell area 
aux.Area = 100 * ones(N*N);
# exponentiate and normalize productivity shocks
aux.Aᵢ = exp.(aux.Aᵢ);
aux.Aᵢ[aux.Country.==1] = aux.Aᵢ[aux.Country.==1] ./ geomean(aux.Aᵢ[aux.Country.==1]);
aux.Aᵢ[aux.Country.==0] = aux.Aᵢ[aux.Country.==0] ./ geomean(aux.Aᵢ[aux.Country.==0]);
[mean(aux.Aᵢ) std(aux.Aᵢ) maximum(aux.Aᵢ) minimum(aux.Aᵢ)];
# back with the geometry
fund = georef(aux, grid);

# plotting figure 1
fig1 = Mke.Figure(size=(600,500))
ax = Mke.Axis(fig1[1, 1], 
            title = "Log productivity",
            titlealign = :left, 
            xlabel = "Longitude", 
            ylabel = "Latitude")
viz!(ax, fund.geometry, color = log.(fund.Aᵢ), colormap = :viridis)
Mke.Colorbar(fig1[1, 2], limits = (minimum(log.(fund.Aᵢ)), maximum(log.(fund.Aᵢ))), 
        colormap = :viridis, flipaxis = true, ticks=-3:1:3, label="Log points", flip_vertical_label=true)
Mke.save("./figures/productivity.png", fig1, px_per_unit = 900/96) #900 dpi

# *********************
# **** Parameters  **** 
# *********************

α = 0.75 ;# Davis & Ortalo-Magne (2011)
σ = 5.0 ;# Simonovska & Waugh (2014)
#ϕ = 0.375 # Based on international trade data
ϕ = 0.33 ;# Based on international trade data
d = dist .^ ϕ ;# trade costs are a constant elasticity function of effective distance by assumption
d[diagind(d)] .= 1 ;# iceberg transport costs are one
τⁱ = ones(N*N,N*N).*2; # # cell trade tax (by assumption)
τⁱ[diagind(τⁱ)] .= 1; # by construction
τⁱ = τⁱ .- 1;
τᵒ = ones(N*N,N*N); # country trade tax (by assumption)
τᵒ[fund.Country.==1,fund.Country.==0] .=2 ;# east to west
τᵒ[fund.Country.==0,fund.Country.==1] .=2; # west to east
τᵒ = τᵒ .- 1;
F = 1.0; # production fixed cost
LL=153889 ;# US civilian labor force 2010 (Statistical Abstract, millions)
LLwest = LL*(sum(fund.Country.==1)/size(values(fund),1));
LLeast = LL*(sum(fund.Country.==0)/size(values(fund),1)) ;

# ****************************************************************************
# **** Open Economy Solve for Endogenous Variables in Initial Equilibrium ****
# ****************************************************************************

println(">>>> Start Wage and Population Convergence <<<<")

@time wᵢ, Lᵢ, πₙᵢ, converged = solveHLwCtyOpen_E(fund,d,τⁱ,τᵒ,N*N);

# prices, rents, and real wage
Pₙ = Hpindex(fund,Lᵢ,wᵢ,πₙᵢ);
rₙ = Hlandprice(fund,Lᵢ,wᵢ); 
V̅ = Hrealw(fund,Lᵢ,πₙᵢ);

# plotting figure 2
eq_plots(Lᵢ,wᵢ,rₙ,Pₙ);

# ************************************************************************
# ***** Counterfactual Eliminating Border Frictions Between Countries ****
# ************************************************************************

τᵒᶜ = ones(N*N,N*N) ; # counterfactual tax regime

println(">>>> Start Wage and Population Convergence <<<<")

@time wᵢᶜ, Lᵢᶜ, πₙᵢᶜ, convergedᶜ =solveHLwCtyOpen_E(fund,d,τⁱ,τᵒᶜ,N*N);

# prices, rents, and real wage
Pₙᶜ = Hpindex(fund,Lᵢᶜ,wᵢᶜ,πₙᵢᶜ);
rₙᶜ = Hlandprice(fund,Lᵢᶜ,wᵢᶜ)  ;
V̅ᶜ = Hrealw(fund,Lᵢᶜ,πₙᵢᶜ);

# exact hat algebra 
L̂ᵢ=Lᵢᶜ./Lᵢ ;
ŵᵢ=wᵢᶜ./wᵢ ;
r̂ₙ=rₙᶜ./rₙ ;
P̂ₙ=Pₙᶜ./Pₙ ;

# exact hat algebra of welfare (welfare gains)
V̂ = Hwelfaregains(πₙᵢ,πₙᵢᶜ,L̂ᵢ) ;
V̂ = round.(V̂, digits=4) ;

# plotting figure 3
eq_plots(L̂ᵢ,ŵᵢ,r̂ₙ,P̂ₙ,suffix="_counterfactual") ;

# **************************************************************************
# ***** Counterfactual Eliminating Border Frictions Between Grid Points ****
# **************************************************************************

τⁱᵪ = ones(N*N,N*N); # counterfactual tax regime

@time wᵢᵪ , Lᵢᵪ, πₙᵢᵪ, convergedᵪ =solveHLwCtyOpen_E(fund,d,τⁱᵪ,τᵒ,N*N);

# prices, rents, and real wage
Pₙᵪ = Hpindex(fund,Lᵢᵪ,wᵢᵪ,πₙᵢᵪ);
rₙᵪ = Hlandprice(fund,Lᵢᵪ,wᵢᵪ)  ;
V̅ᵪ = Hrealw(fund,Lᵢᵪ,πₙᵢᵪ);

# exact hat algebra 
L̂ᵢᵪ=Lᵢᵪ./Lᵢ; 
ŵᵢᵪ=wᵢᵪ./wᵢ ;
r̂ₙᵪ=rₙᵪ./rₙ ;
P̂ₙᵪ=Pₙᵪ./Pₙ;

# exact hat algebra of welfare (welfare gains)
V̂ᵪ = Hwelfaregains(πₙᵢ,πₙᵢᵪ,L̂ᵢᵪ);
V̂ᵪ = round.(V̂ᵪ, digits=4);

# plotting figure 4
eq_plots(L̂ᵢᵪ,ŵᵢᵪ,r̂ₙᵪ,P̂ₙᵪ,suffix="_counterfactual2");

# **********************************************
# ***** Counterfactual Scenarios Comparison ****
# **********************************************
using DelimitedFiles

# Table 1

labs = ["External liberalization" ; "Internal liberalization"]

tab = [round.(unique(V̂[fund.Country.==1].-1)*100,digits=1)[1,1] round.(unique(V̂[fund.Country.==0].-1)*100,digits=1)[1,1];
        round.(unique(V̂ᵪ[fund.Country.==1].-1)*100,digits=1)[1,1] round.(unique(V̂ᵪ[fund.Country.==0].-1)*100,digits=1)[1,1]]

fragment_table(labs,tab,"tab1", path="./tables/")