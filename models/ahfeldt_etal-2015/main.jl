# *********************
# **** Load Files  **** 
# *********************
using  Plots, LaTeXStrings, FixedEffectModels, GeoStats, GeoIO, CSV, DataFrames, Statistics, LinearAlgebra, MAT, Random, StatsBase
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

# *****************
# *** read data ***
# *****************

# ** Random Number ** 
s = MersenneTwister(1)
Random.seed!(s)

# ** set parameters
α=0.80; β=0.75; # Set parameter values to values from the literature                                                                  
ν=κε=0.07; # Set commuting decay to reduced-form estimate

# read .mat data
fileIn = matopen("./data/input/prepdata_big_TD86.mat")
dset = read(fileIn)
close(fileIn)
"
It follows a brief data dictionary:
----------
nobs86 = number of observations
floor86 = rent prices
empwpl86 = workplace employment (population)
emprsd86rw = residential employment (population)
tt86rw = bilateral travel time matrix s.t. rows (i) denote workplaces and columns (j) denote residences
----------
"
S = Int64(dset["nobs86rw"]); Qⱼ = dset["floor86rw"]; Hₘⱼ = dset["empwpl86rw"] ; Hᵣᵢ = dset["emprsd86rw"]; τᵢⱼ = dset["tt86rw"]; # using the paper notation
block_bzk = dset["bzk86rw"] # mapping of Blocks to Bezirkes  

# ** Bezirke (district) wage data ** 
bzkwge = CSV.read("./data/input/wageworker1986.csv", DataFrame; header = false) # read raw data
lwⱼ = log.(bzkwge.Column2) # taking log
lwⱼ = lwⱼ - mean(lwⱼ) # demean wages
Vlwⱼ = var(lwⱼ) # compute variance of log wages, our empirical moment




