module Types

export EndogenousModelParameters, ModelParameters, ExogenousFundamentals, PricesGuess, payroll_aggregator_parameters, EstimationParameters, InverterInputs

using SparseArrays

struct InverterInputs
    Qⱼ::Vector{Float64} # rent prices
    Hₘⱼ::Vector{Float64} # workplace employment
    Hᵣᵢ::Vector{Float64} # residence employment
    τᵢⱼ::Matrix{Float64} # travel time matrix
    Kᵢ::Vector{Float64} # geographical area
end

struct EndogenousModelParameters
    λ::Float64 # Relative importance of production externalities in overall productivity.
    δ::Float64 # Rate of spatial decay (productivity) 
    η::Float64 # Net effect of residence employment density on amenities
    ρ::Float64 # Rate of spatial decay (amenities)
end

struct ModelParameters
    α::Float64 # Share Of consumption
    β::Float64 # Share Of Labor Costs
    μ::Float64 # Share Of Capital In Construction Costs
    κ::Float64 # Size Of Commuting Costs
    ε::Float64 # Fréchet Shape Parameter
end

struct ExogenousFundamentals
    Ãⱼ::Vector{Float64} # Workplace Productivity
    B̃ᵢ::Vector{Float64} # Residential Amenities 
    φᵢ::Vector{Float64} # Density Of Development 
    Kᵢ::Vector{Float64} # Geographical Area
    τᵢⱼ::Matrix{Float64} # Bilateral Travel Time Matrix (i=coming from; j=going to)
end

struct PricesGuess
    Qⱼ0::Vector{Float64} # Rent
    w̃ⱼ0::Vector{Float64} # Wages
    θᵢ0::Vector{Float64} # Share Of Commercial Floorspace
end

function PricesGuess(n_places::Int, pure_res::Vector{Bool}, shared_space::Vector{Bool})
    Q = ones(n_places)
    w = ones(n_places)
    θ = ones(n_places)
    θ[pure_res] .= 0
    θ[shared_space] .= 0.5
    return PricesGuess(Q, w, θ)
end

struct payroll_aggregator_parameters
    S::SparseMatrixCSC{Int, Int}
    Hₘⱼ::Vector{Float64}
    ωⱼ::Vector{Float64}
    Vlwⱼ::Float64
end

struct EstimationParameters
    α::Float64 # Share Of consumption
    ν::Float64 #  ν = εκ is the semi-elasticity of commuting flows wrt travel times
end

end