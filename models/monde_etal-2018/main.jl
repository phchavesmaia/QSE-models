"
MRRH = Monte, Redding, and Rossi-Hansberg
SW = Seidel and Wickerath
"

# *********************
# **** Load Files  **** 
# *********************
using  Plots, FixedEffectModels, GeoStats, GeoIO, CSV, DataFrames, Statistics, LinearAlgebra
import CairoMakie as Mke

function load_dir(dir::String)
    files = readdir(dir)
    for file in files
        if endswith(file, ".jl")
            include(joinpath(dir, file))
        end
    end
end

try 
    cd("D:/Dropbox/learn-julia/qse/models/monde_etal-2018/")
catch
    cd("C:/Users/pedro.maia/Dropbox/learn-julia/qse/models/monde_etal-2018/")
end
load_dir("functions")

# *********************
# **** Parameters  **** 
# *********************

α = 0.7; ε = 4.6; μ = 0.47; δ = 0.38; 
σ = 4; f=1; ν = 0.05; ψ = 0.42; J=401; 

# ********************
# **** Read Data  **** 
# ******************** 

# import csv files 
dataComm = CSV.read("./data/input/commuting_wide.csv", DataFrame; header = true)
dataHous = CSV.read("./data/input/house_prices.csv", DataFrame; header = true)
dataDistance = CSV.read("./data/input/distance_matrix.csv", DataFrame; header = true)
dataArea = CSV.read("./data/input/CountyArea.csv", DataFrame; header = true)
labor_file = CSV.read("./data/input/labor_tidy.csv", DataFrame; header = true)
no_traffic = CSV.read("./data/input/roundtrip_time_base.csv", DataFrame; header = true)

# commuting matrices
comMat = Matrix(dataComm[:,2:end])
comMat = comMat'
diffe = sum(comMat, dims=2) .- sum(comMat, dims=1)'
L = sum(comMat)
λₙᵢ = comMat ./ L
λₙᵢi = λₙᵢ ./ sum(λₙᵢ, dims = 2)

# wage data
wₙ = labor_file.median_income_workplace
wₙ = wₙ ./ mean(wₙ) # normalization
v̄ₙ = λₙᵢi * wₙ # expected residential income (eq. 14 in MRRH or eq. 6 in SW) -- this is in the 'model inversion' section

# demographics data
Lₙ = sum(λₙᵢ, dims=1)' .* L # people who work at n -- this is in the 'model inversion' section
Lₙ = Lₙ ./ mean(Lₙ) # normalization
Rₙ = sum(λₙᵢ, dims=2) .* L # people who live in n -- this is in the 'model inversion' section
Rₙ = Rₙ ./ mean(Rₙ) # normalization

# housing prices data
Pₕₙ = dataHous.rentindex
lPₕₙ = log.(Pₕₙ)

# distance data and computing trade costs from distances
distₙᵢ =  Matrix(dataDistance[:,2:end])
distₙᵢ = distₙᵢ ./ minimum(distₙᵢ)
dₙᵢ = distₙᵢ .^ ψ # nonumbered equation in both papers, but its on the "model inversion" section (fund. produc.)

# geographic area
Areaₙ = dataArea.Area

# congestion (traffic) data
baseline = Matrix(no_traffic[:,2:end])
baseline = baseline'

# *************************************************************
# **** Productivity and Trade Shares from Model Inversion  **** 
# *************************************************************

@time Aᵢ, πₙᵢ, Pₙ  = solveProductTrade(Lₙ, Rₙ, wₙ, v̄ₙ, dₙᵢ)
println("<<<<<<<<<<<<<<< Data compilation completed >>>>>>>>>>>>>>>")

# *******************************
# **** Descriptive Analysis  **** 
# *******************************

descriptive_analysis(Lₙ, Areaₙ, Rₙ, distₙᵢ, comMat, Aᵢ)

# **************************
# **** Counterfactuals  **** 
# **************************
#= 
    Since this counterfatual analysis is not explained in neither research papers,
    it makes sense to have a brief discussion of what I am trying to implement. The
    goal of this counterfactual exercise is to assess the GE effects of hypothetical
    border that follows the border between formerly separated West Germany and East 
    Germany.
=#
# ----------------- Border Data Prep --------------------

# Create a dummy for counties in the East
east = zeros(J,1)
east[325:end] .= 1 # Set elements from 325 to the end to 1 

# Create the condition matrix using outer logical OR between East-West and West-East
conditionMatrix = ((east.==1) * (east.==0)') .+ ((east.==0) * (east.==1)')

# Import border distance
dataBoderDist = CSV.read("./data/input/CountyBorderDist.csv", DataFrame; header = true)
BorderDist_n = dataBoderDist.BorderDist .+ 10 # Add a small distance to improve visibility in scatter plot
BorderDist_n[east.==0] = -abs.(BorderDist_n[east.==0]) # value 0 in the center of germany from west (negative) to east (positive)

# -------------- Counterfactual Estimation --------------

