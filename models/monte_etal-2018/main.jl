"
MRRH = Monte, Redding, and Rossi-Hansberg
SW = Seidel and Wickerath
ARSW = Afhdelt et. al
"

# *********************
# **** Load Files  **** 
# *********************
using  Plots, LaTeXStrings, FixedEffectModels, GeoStats, GeoIO, CSV, DataFrames, Statistics, LinearAlgebra
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
    cd("D:/Dropbox/learn-julia/qse/models/monte_etal-2018/")
catch
    cd("C:/Users/pedro.maia/Dropbox/learn-julia/qse/models/monte_etal-2018/")
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
dataComm = CSV.read("./data/input/commuting_wide.csv", DataFrame; header = true);
dataHous = CSV.read("./data/input/house_prices.csv", DataFrame; header = true);
dataDistance = CSV.read("./data/input/distance_matrix.csv", DataFrame; header = true);
dataArea = CSV.read("./data/input/CountyArea.csv", DataFrame; header = true);
labor_file = CSV.read("./data/input/labor_tidy.csv", DataFrame; header = true);
no_traffic = CSV.read("./data/input/roundtrip_time_base.csv", DataFrame; header = true);

# commuting matrices (Assumes you actually have the commuting flows matrix)
comMat = Matrix(dataComm[:,2:end]);
comMat = comMat';
diffe = sum(comMat, dims=2) .- sum(comMat, dims=1)';
L = sum(comMat);
λₙᵢ = comMat ./ L;
λₙᵢn = λₙᵢ ./ sum(λₙᵢ, dims = 2);

# wage data
wₙ = labor_file.median_income_workplace;
wₙ = wₙ ./ mean(wₙ); # normalization
v̄ₙ = λₙᵢn * wₙ ;# expected residential income (eq. 14 in MRRH or eq. 6 in SW)

# demographics data
Lₙ = sum(λₙᵢ, dims=1)' .* L; # people who work at n -- this is in the 'model inversion' section
Lₙ = Lₙ ./ mean(Lₙ); # normalization
Rₙ = sum(λₙᵢ, dims=2) .* L ; # people who live in n -- this is in the 'model inversion' section
Rₙ = Rₙ ./ mean(Rₙ); # normalization

# housing prices data
Pₕₙ = dataHous.rentindex;
lPₕₙ = log.(Pₕₙ);

# distance data and computing trade costs from distances
distₙᵢ =  Matrix(dataDistance[:,2:end]);
distₙᵢ = distₙᵢ ./ minimum(distₙᵢ);
dₙᵢ = distₙᵢ .^ ψ ;# nonumbered equation in both papers, but it is on the "model inversion" section (fund. produc.)

# geographic area
Areaₙ = dataArea.Area;

# congestion (traffic) data
baseline = Matrix(no_traffic[:,2:end]);
baseline = baseline';


@time Tᵢ, λₙᵢn_model, λₙᵢ_model, Lₙᵢ = get_bi_lamba(wₙ, baseline, Rₙ, Lₙ , L);  # The function gets λₙᵢn and λₙᵢ if you don't have it in your data. 
                                                                             # If you do not have wₙ, let wₙ=1 and interpret transformed wages ωₙ = Tᵢ (close to eq. 26 of ARSW). 
                                                                             # Then, we have that model-consistent wages can then be recovered as w̃ᵢ = Tᵢ .^ (1/ε) 

# Let's compare the estimated λₙᵢ_model with the real λₙᵢ...
df = DataFrame()
df[!,"real"] = vec(reshape(λₙᵢ,:,1))
df[!,"model"] = vec(reshape(λₙᵢ_model,:,1))
model = reg(df, @formula(real ~ model))
scatter(reshape(λₙᵢ,:,1), reshape(λₙᵢ_model,:,1), 
        legend = false, xlabel = "Real λₙᵢ", ylabel="Model consistent λₙᵢ",
        title="Model-consistent λₙᵢ fit", smooth=:true)
annotate!(
    0.01,
    0.035,
    latexstring("R^2 = $(round(r2(model), digits = 2))")
    )
savefig("./figures/model_consistent_lambda_fit.png");

# ************************************************************************
# **** Productivity, Amenities and Trade Shares from Model Inversion  **** 
# ************************************************************************

@time Aᵢ, πₙᵢ, Pₙ  = solveProductTrade(Lₙ, Rₙ, wₙ, v̄ₙ, dₙᵢ);
println("<<<<<<<<<<<<<<< Data compilation completed >>>>>>>>>>>>>>>")

# *******************************
# **** Descriptive Analysis  **** 
# *******************************

descriptive_analysis(Lₙ, Areaₙ, Rₙ, distₙᵢ, comMat, Aᵢ, baseline, λₙᵢ);

# **************************
# **** Counterfactuals  **** 
# **************************
"
    Since this counterfatual analysis is not explained in neither research papers,
    it makes sense to have a brief discussion of what I am trying to implement. The
    goal of this counterfactual exercise is to assess the GE effects of hypothetical
    border that follows the border between formerly separated West Germany and East 
    Germany.
"
# ----------------- Border Data Prep --------------------

# Create a dummy for counties in the East
east = zeros(J,1);
east[325:end] .= 1 ;# Set elements from 325 to the end to 1 

# Create the condition matrix using outer logical OR between East-West and West-East
conditionMatrix = ((east.==1) * (east.==0)') .+ ((east.==0) * (east.==1)');

# Import border distance
dataBoderDist = CSV.read("./data/input/CountyBorderDist.csv", DataFrame; header = true);
BorderDist_n = dataBoderDist.BorderDist .+ 10; # Add a small distance to improve visibility in scatter plot
BorderDist_n[east.==0] = -abs.(BorderDist_n[east.==0]) ;# value 0 at the germany division point from west (negative) to east (positive)

# -------------- Creating a *trade* border between East and West Germany --------------
Âᵢ=ones(size(wₙ)); B̂ₙᵢ=ones(size(πₙᵢ)); κ̂ₙᵢ=ones(size(πₙᵢ));
# creating the tax  
d̂ₙᵢ = ones(size(dₙᵢ));
d̂ₙᵢ[conditionMatrix.==1] .= 1000 ;# Set border cost in trade to large value

# estimating counterfactual variables 
@time ŵₙ, v̄̂ₙ, P̂ₕₙ, π̂ₙᵢ, λ̂ₙᵢ, P̂ₙ, R̂ₙ, L̂ₙ, Ū̂ = counterFacts(wₙ,v̄ₙ,Lₙ,Rₙ,πₙᵢ,λₙᵢ, d̂ₙᵢ=d̂ₙᵢ);

# Map findings
mapit("./data/shapes/VG250_KRS_clean_final.shp",log.(P̂ₕₙ),"Relative change in house price; trade border", path_to="./figures/MAP_COUNT_BORDER_TRADE_Qchange.png") ;
mapit("./data/shapes/VG250_KRS_clean_final.shp",log.(P̂ₙ),"Relative change in tradable goods price; trade border", path_to="./figures/MAP_COUNT_BORDER_TRADE_Pchange.png") ;
mapit("./data/shapes/VG250_KRS_clean_final.shp",log.(R̂ₙ),"Relative change in population; trade border", path_to="./figures/MAP_COUNT_BORDER_TRADE_Rchange.png") ;

"
Essentially, the two sides of Germany simply stopped trading with 
each other, hence prices go downs due to the loss of potential market,
i.e., loss of (feasible) demand for the tradable good. This implies
less firms, leading to less jobs (lower R̂ₙ in equilibrium) and lower 
wages, which shift the non-tradable good (housing) prices down.
"

# -------------- Creating a *commuting* border between East and West Germany --------------

# destroying every commuting road between east and west germany 
κ̂ₙᵢ = ones(size(dₙᵢ));
κ̂ₙᵢ[conditionMatrix.==1] .= 1000;

# estimating counterfactual variables 
@time ŵₙ, v̄̂ₙ, P̂ₕₙ, π̂ₙᵢ, λ̂ₙᵢ, P̂ₙ, R̂ₙ, L̂ₙ, Ū̂ = counterFacts(wₙ,v̄ₙ,Lₙ,Rₙ,πₙᵢ,λₙᵢ, κ̂ₙᵢ=κ̂ₙᵢ);

# Map findings
mapit("./data/shapes/VG250_KRS_clean_final.shp",log.(P̂ₕₙ),"Relative change in house price; trade border", path_to="./figures/MAP_COUNT_BORDER_COMM_Qchange.png"); 
mapit("./data/shapes/VG250_KRS_clean_final.shp",log.(P̂ₙ),"Relative change in tradable goods price; trade border", path_to="./figures/MAP_COUNT_BORDER_COMM_Pchange.png") ;
mapit("./data/shapes/VG250_KRS_clean_final.shp",log.(R̂ₙ),"Relative change in population; trade border", path_to="./figures/MAP_COUNT_BORDER_COMM_Rchange.png") ;

# -------------- Creating a *commuting* and *trade* border border between East and West Germany --------------

# estimating counterfactual variables 
@time ŵₙ, v̄̂ₙ, P̂ₕₙ, π̂ₙᵢ, λ̂ₙᵢ, P̂ₙ, R̂ₙ, L̂ₙ, Ū̂ = counterFacts(wₙ,v̄ₙ,Lₙ,Rₙ,πₙᵢ,λₙᵢ, d̂ₙᵢ=d̂ₙᵢ, κ̂ₙᵢ=κ̂ₙᵢ);

# Map findings
mapit("./data/shapes/VG250_KRS_clean_final.shp",log.(P̂ₕₙ),"Relative change in house price; trade border", path_to="./figures/MAP_COUNT_BORDER_COMMTRADE_Qchange.png") ;
mapit("./data/shapes/VG250_KRS_clean_final.shp",log.(P̂ₙ),"Relative change in tradable goods price; trade border", path_to="./figures/MAP_COUNT_BORDER_COMMTRADE_Pchange.png") ;
mapit("./data/shapes/VG250_KRS_clean_final.shp",log.(R̂ₙ),"Relative change in population; trade border", path_to="./figures/MAP_COUNT_BORDER_COMMTRADE_Rchange.png") ;

# Scatter plot of goods market impact 
scatter(BorderDist_n, log.(P̂ₕₙ), 
        xlabel="Border Distance (km)", 
        ylabel="Log change in prices", 
        label = "Housing", 
        title="Effects of introducing a domestic border (prices)", 
        grid=true, 
        legend=:bottomleft)
scatter!(BorderDist_n, log.(P̂ₙ), label = "Tradable goods", grid=true)
savefig("./figures/scatter_COUNT_BORDER_COMMTRADE_PriceChanges.png");

# Scatter plot of labour market impact 
scatter(BorderDist_n, log.(ŵₙ), 
        xlabel="Border Distance (km)", 
        ylabel="Log change", 
        label = "Wages", 
        title="Effects of introducing a domestic border (labour)", 
        grid=true, 
        legend=:bottomleft)
scatter!(BorderDist_n, log.(L̂ₙ), label = "Employment", grid=true)
savefig("./figures/scatter_COUNT_BORDER_COMMTRADE_LabourChanges.png");


# -------------- Now repeat with half the trade cost --------------
ψ = 0.21; dₙᵢ = distₙᵢ .^ ψ; # lowering trade costs

# Productivity and trade shares from model inversion
@time Aᵢ, πₙᵢ, Pₙ  = solveProductTrade(Lₙ, Rₙ, wₙ, v̄ₙ, dₙᵢ);

# estimating counterfactual variables 
@time ŵₙ, v̄̂ₙ, P̂ₕₙ, π̂ₙᵢ, λ̂ₙᵢ, P̂ₙ, R̂ₙ, L̂ₙ, Ū̂ = counterFacts(wₙ,v̄ₙ,Lₙ,Rₙ,πₙᵢ,λₙᵢ, d̂ₙᵢ=d̂ₙᵢ, κ̂ₙᵢ=κ̂ₙᵢ);

# Map findings
mapit("./data/shapes/VG250_KRS_clean_final.shp",log.(P̂ₕₙ),"Relative change in house price; trade border", path_to="./figures/MAP_COUNT_BORDER_COMMTRADE_Qchange_lowTradeCost.png"); 
mapit("./data/shapes/VG250_KRS_clean_final.shp",log.(P̂ₙ),"Relative change in tradable goods price; trade border", path_to="./figures/MAP_COUNT_BORDER_COMMTRADE_Pchange_lowTradeCost.png");
mapit("./data/shapes/VG250_KRS_clean_final.shp",log.(R̂ₙ),"Relative change in population; trade border", path_to="./figures/MAP_COUNT_BORDER_COMMTRADE_Rchange_lowTradeCost.png");

# Scatter plot of goods market impact 
scatter(BorderDist_n, log.(P̂ₕₙ), 
        xlabel="Border Distance (km)", 
        ylabel="Log change in prices", 
        label = "Housing", 
        title="Effects of introducing a domestic border (prices)", 
        grid=true, 
        legend=:bottomleft)
scatter!(BorderDist_n, log.(P̂ₙ), label = "Tradable goods", grid=true)
savefig("./figures/scatter_COUNT_BORDER_COMMTRADE_PriceChanges_lowTradeCost.png");

# Scatter plot of labour market impact 
scatter(BorderDist_n, log.(ŵₙ), 
        xlabel="Border Distance (km)", 
        ylabel="Log change", 
        label = "Wages", 
        title="Effects of introducing a domestic border (labour)", 
        grid=true, 
        legend=:bottomleft)
scatter!(BorderDist_n, log.(L̂ₙ), label = "Employment", grid=true)
savefig("./figures/scatter_COUNT_BORDER_COMMTRADE_LabourChanges_lowTradeCost.png");


