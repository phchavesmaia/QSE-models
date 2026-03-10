# *********************
# **** Load Files  **** 
# *********************
using  Plots, LaTeXStrings, FixedEffectModels, GeoStats, GeoIO, CSV, DataFrames, Statistics, LinearAlgebra, MAT, Random, StatsBase, Optimization, OptimizationNLopt, BenchmarkTools, SparseArrays
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
    cd("/home/phchavesmaia/Dropbox/learn-julia/qse/models/ahfeldt_etal-2015")
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

# ****************************************************************
# *** Estimating the frechet-shape parameter ε using 1986 data ***
# ****************************************************************
"
    This section of the code aims to compute the value of the
    structural parameter ε. We validade this estimate by 
    analyzing wheter it would lead to model estimates of 
    workforce allocation consistent with 'real world' data.
"
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
S = Int64(dset["nobs86rw"]); Qⱼ = dset["floor86rw"]; Hₘⱼ = dset["empwpl86rw"] ; Hᵣᵢ = dset["emprsd86rw"]; τᵢⱼ = dset["tt86rw"]'; # using the paper notation
block_bzk = dset["bzk86rw"]; # map from blocks to districts 
dset = nothing; GC.gc() # free memory

bzkwge = CSV.read("./data/input/wageworker1986.csv", DataFrame; header = false); # Bezirke (district) raw wage data
lwⱼ = log.(bzkwge.Column2); # taking log
lwⱼ = lwⱼ .- mean(lwⱼ); # demean wages
Vlwⱼ = var(lwⱼ); # compute variance of log wages, our empirical moment

# computing ω and ε using 86 data
ε, Ĥₘⱼ, ωⱼ = get_ε(Vlwⱼ,Hₘⱼ,Hᵣᵢ,τᵢⱼ,Qⱼ, su=block_bzk); 
ε = round(ε, digits=2); # rounded for consistency with replication package
κ = round(ν/ε,digits=6); # setting commuting decay to reduced-form estimate; rounded for consistency with replication package

# Defining a sanity check function
function snty_check(v1,v2; tol=6)
    df = DataFrame()
    df[!,"model"] = vec(v1)
    df[!,"data"] = vec(v2)
    model = reg(df,@formula(data ~ model));
    slope = round(coef(model)[2],digits=tol);
    return slope
end

# 'validating' estimates
println("The slope of the model/real workplace population data is: $(snty_check(Hₘⱼ,Ĥₘⱼ))")

# *********************************************************************
# *** Calibration of exogenous fundamentals (model inversion, 2006) ***
# *********************************************************************
"
    As explained in the ARSW Codebook, the functions implemented in this section
    aim to recover the structural fundamentals {Ãᵢ,B̃ᵢ,φᵢ}, which incorporate 
    {Tᵢ,Eᵢ,Aᵢ,Bᵢ,φᵢ,Kᵢ,ξᵢ}, by inverting the model so that we can use endogenous 
    variables to recover the structural parameters that rationalize them.
"

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
Qⱼ = dset["floor06"]; Hₘⱼ = dset["empwpl06"] ; Hᵣᵢ = dset["emprsd06"]; τᵢⱼ = dset["tt06"]; Kᵢ = dset["area06"];
block_bzk06 = dset["bzk06"]; 
dset = nothing; GC.gc() # free memory

# computing the structural fundamentals of the model SEQUENTIALLY
Ãⱼ, B̃ᵢ, w̃ⱼ, πᵢⱼ, Tw̃ᵢ, φᵢ, Lᵢ, θᵢ, H̃ₘⱼ, H̃ᵣᵢ, CMA = cal_model_seq(Qⱼ,Hₘⱼ,Hᵣᵢ,τᵢⱼ,Kᵢ); 

# computing the structural fundamentals of the model SIMULTANEOUSLY
Ãⱼsim, B̃ᵢsim, w̃ⱼsim, πᵢⱼsim, Tw̃ᵢsim, φᵢsim, Lᵢsim, θᵢsim, H̃ₘⱼsim, H̃ᵣᵢsim, CMAsim = cal_model_sim(Qⱼ,Hₘⱼ,Hᵣᵢ,τᵢⱼ,Kᵢ); 

# sanity check that both algorithms derive the same results
snty_check_algo = [
    snty_check(Ãⱼ, Ãⱼsim), 
    snty_check(B̃ᵢ, B̃ᵢsim),
    snty_check(w̃ⱼ, w̃ⱼsim),
    snty_check(πᵢⱼ, πᵢⱼsim),
    snty_check(Tw̃ᵢ, Tw̃ᵢsim), 
    snty_check(φᵢ, φᵢsim), 
    snty_check(Lᵢ, Lᵢsim), 
    snty_check(θᵢ, θᵢsim), 
    snty_check(H̃ₘⱼ, H̃ₘⱼsim), 
    snty_check(H̃ᵣᵢ, H̃ᵣᵢsim), 
    snty_check(CMA, CMAsim, tol=5)];
println("Are the sequential and simultaneous algorithms agreeing? $(sum(snty_check_algo)==length(snty_check_algo))")

"--- In benchmark tests (@btime from a cold start) comparing the SEQUENTIAL with the SIMULTANEOUS approaches, I have found that:
 ---     1. SEQUENTIAL takes: 45.691 s (54603962 allocations: 90.55 GiB)
 ---     2. SIMULTANEOUS takes: 424.806 s (4393 allocations: 6.49 GiB)
 --- Since: 
 ---     1) both algorithms return virtually identical results; and 
 ---     2) Julia's garbage collects is quite impressive 
 --- I would highly recommend opting for the SEQUENTIAL approach.
"

# indirectly 'validating' bilateral commuting probabilities calibration (it should be 1!)
println("The slope of the model-real workplace population data is: $(snty_check(Hₘⱼ,H̃ₘⱼ))")

# plotting some maps
mapit("./data/shapefile/Berlin4matlab1.shp",CMA,"Commuter market access", label_legend="", path_to="./figures/cma06.png");
mapit("./data/shapefile/Berlin4matlab1.shp",Ãⱼ,"Productivity", label_legend="", path_to="./figures/productivity06.png");
mapit("./data/shapefile/Berlin4matlab1.shp",B̃ᵢ,"Amenities", label_legend="", path_to="./figures/amenities06.png");
mapit("./data/shapefile/Berlin4matlab1.shp",φᵢ,"Density of Development", label_legend="", path_to="./figures/density06.png");

# *********************************************************
# *** Solving 2006 equilibrium (exogenous fundamentals) *** 
# *********************************************************
"
    This section of the code will solve the model by inputing
    the structural parameters and exogenous fundamentals in
    order to recover the equilibrium endogenous variables.
    Notably, we must use the {Ãᵢ,B̃ᵢ,φᵢ} estimates recovered
    in the previous section. 
"
# enunciate model parameters
params = (α, β, κ, ε, μ);

# enunciate exogenous fundamentals of the model
exo_fund = (Ãⱼ, B̃ᵢ, φᵢ, Kᵢ, τᵢⱼ); 

# enunciate guesses at equilibrium prices
prices_guess = (Qⱼ, w̃ⱼ, θᵢ); 

# solve the equilibrium using data/model-consistent initial guesses
w̃ⱼeq, θᵢeq, Qⱼeq, πᵢⱼeq, H̃eq = solve_equilibrium(params, exo_fund, prices_guess = prices_guess);

# validating the equilibrium variables with real data
snty_check_eq = [
    snty_check(w̃ⱼeq,w̃ⱼ),
    snty_check(θᵢeq,θᵢ),
    snty_check(Qⱼeq,Qⱼ),
    snty_check(πᵢⱼeq,πᵢⱼ),
    Int(round(H̃eq,digits=0) == round(sum(Hₘⱼ),digits=0))
];
println("Does this equilibrium match the data? $(sum(snty_check_eq)==length(snty_check_eq))")

# solve the equilibrium blindly (it takes a long time!)
w̃ⱼeqb, θᵢeqb, Qⱼeqb, πᵢⱼeqb, H̃eqb = solve_equilibrium(params, exo_fund)


