# *********************
# **** Load Files  **** 
# *********************
using  LaTeXStrings, FixedEffectModels, CSV, DataFrames, Statistics, Random, BenchmarkTools, MAT 

try 
    cd("/home/phchavesmaia/Dropbox/learn-julia/qse/models/ahfeldt_etal-2015")
catch
    cd("C:/Users/pedro.maia/Dropbox/learn-julia/qse/models/ahfeldt_etal-2015/")
end
include("./modules/types.jl")
include("./modules/matlab_helpers.jl")
include("./modules/frechet_estimation.jl")
include("./modules/model_inverters.jl")
include("./modules/model_solvers.jl")

# **************************
# *** Setting parameters ***
# **************************

# define model parameters
module Parameters
    # Set parameter values to values from the literature
    const α=0.80; 
    const β=0.75; 
    const μ=0.75; 
    # Set commuting decay to reduced-form estimate (see ARSW table 3)
    const ν=κε=0.07; 
    # export parameters
    export α, β, μ, ν, κε
end 

using .Parameters, .Types, .MatlabHelpers, .FrechetEstimation, .ModelInverters, .ModelSolver

# Random Number 
s = MersenneTwister(1);
Random.seed!(s);

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
Qⱼ, Hₘⱼ, Hᵣᵢ, τᵢⱼ, Kᵢ, block_bzk = read_mat("86");

bzkwge = CSV.read("./data/input/wageworker1986.csv", DataFrame; header = false); # Bezirke (district) raw wage data
lwⱼ = @. log(bzkwge.Column2); # taking log
lwⱼ = @. lwⱼ - $mean(lwⱼ); # demean wages
Vlwⱼ = var(lwⱼ); # compute variance of log wages, our empirical moment

# computing ω and ε using 86 data
params = EstimationParameters(α,ν);
ε⁼, Ĥₘⱼ, ωⱼ = get_ε(Vlwⱼ,Hₘⱼ,Hᵣᵢ,τᵢⱼ,Qⱼ,params, su=block_bzk); 
const ε = round(ε⁼, digits=2); # rounded for consistency with replication package
const κ = round(ν/ε,digits=6); # setting commuting decay to reduced-form estimate; rounded for consistency with replication package

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

GC.gc() # garbage collector (free memory)

"
    As explained in the ARSW Codebook, the functions implemented in this section
    aim to recover the structural fundamentals {Ãᵢ,B̃ᵢ,φᵢ}, which incorporate 
    {Tᵢ,Eᵢ,Aᵢ,Bᵢ,φᵢ,Kᵢ,ξᵢ}, by inverting the model so that we can use endogenous 
    variables to recover the structural parameters that rationalize them.
"

# read 2006 data 
Qⱼ, Hₘⱼ, Hᵣᵢ, τᵢⱼ, Kᵢ, block_bzk06 = read_mat("06");

# enunciate model parameters
params = ModelParameters(α, β, κ, ε, μ);
inputs = InverterInputs(Qⱼ,Hₘⱼ,Hᵣᵢ,τᵢⱼ,Kᵢ);

# computing the structural fundamentals of the model SEQUENTIALLY
Ãⱼ, B̃ᵢ, w̃ⱼ, πᵢⱼ, Tw̃ᵢ, φᵢ, Lᵢ, θᵢ, H̃ₘⱼ, H̃ᵣᵢ, CMA = invert_model(params, inputs); 

# computing the structural fundamentals of the model SIMULTANEOUSLY
Ãⱼsim, B̃ᵢsim, w̃ⱼsim, πᵢⱼsim, Tw̃ᵢsim, φᵢsim, Lᵢsim, θᵢsim, H̃ₘⱼsim, H̃ᵣᵢsim, CMAsim = invert_model(params, inputs, method="simultaneous", stop_rule="matlab"); 

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

GC.gc() # garbage collector (free memory)

"
    This section of the code will solve the model by inputing
    the structural parameters and exogenous fundamentals in
    order to recover the equilibrium endogenous variables.
    Notably, we must use the {Ãᵢ,B̃ᵢ,φᵢ} estimates recovered
    in the previous section. 
"

# enunciate exogenous fundamentals of the model
exo_fund = ExogenousFundamentals(Ãⱼ, B̃ᵢ, φᵢ, Kᵢ, τᵢⱼ); 

# enunciate guesses at equilibrium prices
prices_guess = PricesGuess(Qⱼ, w̃ⱼ, θᵢ); 

# solve the CLOSED-CITY equilibrium using data/model-consistent initial guesses
H = sum(Hᵣᵢ);
Qⱼeq, w̃ⱼeq, θᵢeq, πᵢⱼeq, Ūeq = solve_equilibrium(params, exo_fund, H, closed_city = true, prices_guess = prices_guess);

# validating the CLOSED-CITY equilibrium variables with real data
snty_check_eq_closed = [
    snty_check(w̃ⱼeq,w̃ⱼ,tol=3),
    snty_check(θᵢeq,θᵢ,tol=3),
    snty_check(Qⱼeq,Qⱼ,tol=3),
    snty_check(πᵢⱼeq,πᵢⱼ,tol=3)];
println("Does this CLOSED-CITY equilibrium match the data? $(sum(snty_check_eq_closed)==length(snty_check_eq_closed))")

# solve THE OPEN-CITY equilibrium using data/model-consistent initial guesses
Qⱼeq_open, w̃ⱼeq_open, θᵢeq_open, πᵢⱼeq_open, H̃eq_open = solve_equilibrium(params, exo_fund, Ūeq, closed_city = false, prices_guess = prices_guess);

# validating the OPEN-CITY equilibrium variables with real data
snty_check_eq_open = [
    snty_check(w̃ⱼ,w̃ⱼeq_open,tol=3),
    snty_check(θᵢ,θᵢeq_open,tol=3),
    snty_check(Qⱼ,Qⱼeq_open,tol=3),
    snty_check(πᵢⱼ,πᵢⱼeq_open,tol=3),
    Int(round(H/H̃eq_open,digits=3))];
    
println("Does this OPEN-CITY equilibrium match the data? $(sum(snty_check_eq_open)==length(snty_check_eq_open))")

# solve the CLOSED-CITY equilibrium blindly (no initial guesses at prices)
Qⱼeqb, w̃ⱼeqb, θᵢeqb, πᵢⱼeqb, Ūeqb = solve_equilibrium(params, exo_fund, H, closed_city = true);

# validaring if initial guesses lead to different results...
snty_check_guess = [
    snty_check(w̃ⱼeq,w̃ⱼeqb,tol=3),
    snty_check(θᵢeq,θᵢeqb,tol=3),
    snty_check(Qⱼeq,Qⱼeqb,tol=3),
    snty_check(πᵢⱼeq,πᵢⱼeqb,tol=3),
    Int(round(Ūeq/Ūeqb,digits=3))];
println("Are the CLOSED-CITY results robust to different initial guesses (is the eq. perfectly identified)? $(sum(snty_check_eq)==length(snty_check_eq))")

# *****************************************************************
# *** Counterfactuals 2006 equilibrium (exogenous fundamentals) *** 
# *****************************************************************

GC.gc() # garbage collector (free memory)

# --- What happens if we ban cars in the entire city? --- #
fileIn = matopen("./data/input/ttpublic_2006_ren.mat");
dset = read(fileIn); close(fileIn);
τᵢⱼpub = Matrix(dset["ttpub06"]); dset = nothing # read counterfactual bilateral travel time matrix

# enunciate altered exogenous fundamentals of the model
exo_fund_ctf = ExogenousFundamentals(Ãⱼ, B̃ᵢ, φᵢ, Kᵢ, τᵢⱼpub); 

# estimate alternative CLOSED-CITY equilibrium
Qⱼpub, w̃ⱼpub, θᵢpub, πᵢⱼpub, Ūpub = solve_equilibrium(params, exo_fund_ctf, H, closed_city = true, prices_guess = prices_guess, tol_digits=2, damp_fact=.4);
println("The welfare change from banning cars would be of: $(round(100*(Ūpub-Ūeq)/Ūeq,digits=2))%")

# estimate alternative OPEN-CITY equilibrium
Qⱼpub_open, w̃ⱼpub_open, θᵢpub_open, πᵢⱼpub_open, H̃pub = solve_equilibrium(params, exo_fund_ctf, Ūeq, closed_city = false, prices_guess = prices_guess, tol_digits=2);
println("The population change from banning cars would be of: $(round(100*(H̃pub/H-1),digits=2))%")

