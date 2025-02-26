# *********************
# **** Load Files  **** 
# *********************
using  Plots, LaTeXStrings, FixedEffectModels, GeoStats, GeoIO, CSV, DataFrames, Statistics, LinearAlgebra, MAT, Random, StatsBase, NLopt, BenchmarkTools
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
    cd("D:/Dropbox/learn-julia/qse/models/ahfeldt_etal-2015/")
catch
    cd("C:/Users/pedro.maia/Dropbox/learn-julia/qse/models/ahfeldt_etal-2015/")
end
load_dir("functions")

# **************************
# *** Setting parameters ***
# **************************

# Random Number 
s = MersenneTwister(1)
Random.seed!(s)

# set parameters
α=0.80; β=0.75; μ=0.75; # Set parameter values to values from the literature                                                                  
ν=κε=0.07; # Set commuting decay to reduced-form estimate (see ARSW table 3)

# *****************************************************
# *** computing ω and the frechet-shape parameter ε ***
# *****************************************************

# read 1986 data
fileIn = matopen("./data/input/prepdata_big_TD86.mat");
dset = read(fileIn);
close(fileIn);
"
It follows a brief data dictionary:
----------
nobs86 = number of observations
floor86 = rent prices
empwpl86 = workplace employment (population)
emprsd86rw = residential employment (population)
tt86rw = bilateral travel time matrix s.t. rows (i) denote workplaces and columns (j) denote residences
bzk86rw = mapping of Blocks to Bezirkes
----------
"
S = Int64(dset["nobs86rw"]); Qⱼ = dset["floor86rw"]; Hₘⱼ = dset["empwpl86rw"] ; Hᵣᵢ = dset["emprsd86rw"]; τᵢⱼ = dset["tt86rw"]; # using the paper notation
block_bzk = dset["bzk86rw"]; # map from blocks to districts 
bzkwge = CSV.read("./data/input/wageworker1986.csv", DataFrame; header = false); # Bezirke (district) raw wage data
lwⱼ = log.(bzkwge.Column2); # taking log
lwⱼ = lwⱼ .- mean(lwⱼ); # demean wages
Vlwⱼ = var(lwⱼ); # compute variance of log wages, our empirical moment

# computing ω and ε using 86 data
ε, Ĥₘⱼ = get_ε(Vlwⱼ,Hₘⱼ,Hᵣᵢ,τᵢⱼ,Qⱼ); 
κ = ν/ε; # setting commuting decay to reduced-form estimate

# 'validating' estimates
df = DataFrame()
df[!,"model"] = vec(Ĥₘⱼ)
df[!,"data"] = vec(Hₘⱼ)
reg(df,@formula(data ~ model))

# ********************************************************************************
# *** Calibration of exogenous fundamentals (solving the equilibrium for 2006) ***
# ********************************************************************************

# read 2006 data
fileIn = matopen("./data/input/prepdata_big_TD.mat");
dset = read(fileIn);
close(fileIn);
"
It follows a brief data dictionary:
----------
nobs06 = number of observations
floor06 = rent prices
empwpl06 = workplace employment (population)
emprsd06 = residential employment (population)
tt06 = bilateral travel time matrix s.t. rows (i) denote workplaces and columns (j) denote residences
bzk06 = mapping of Blocks to Bezirkes
area06 = geographical area 
----------
"
S = Int64(dset["nobs06"]); Qⱼ = dset["floor06"]; Hₘⱼ = dset["empwpl06"] ; Hᵣᵢ = dset["emprsd06"]; τᵢⱼ = dset["tt06"]; Kᵢ = dset["area06"];
block_bzk06 = dset["bzk06"];

# computing structural fundamentals of the model
Ãⱼ, B̃ᵢ, w̃ⱼ, πᵢⱼ, Tw̃ᵢ, ϕᵢ, Lᵢᴰ, θᵢ, H̃ₘⱼ, H̃ᵣᵢ, CMA = cal_model(Qⱼ,Hₘⱼ,Hᵣᵢ,τᵢⱼ,Kᵢ); 