module ModelSolver

using SpecialFunctions 
using ..Types: ModelParameters, EndogenousModelParameters, ExogenousFundamentals, PricesGuess

export solve_equilibrium, detangle_agglomeration

"""
This function assumes that you have predefined the parameters
{α, β, κ, ε, μ}. If you are in the endogenous agglomeration
forces case, you will also need {λ, δ, η, ρ}. It then solves 
for the general equilibrium of the model given the structural 
exogenous fundamentals:

    1. Ãⱼ (or aⱼ) = workplace exogenous productivity;
    2. B̃ᵢ (or bᵢ) = residential exogenous amenities;
    3. φᵢ = density of development;
    4. Kᵢ = geographical area; and
    5. τᵢⱼ = bilateral travel time matrix.
    
The output of this function is the set of endogenous equilibrium
variables {w̃ⱼ, θᵢ, Qⱼ, πᵢⱼ} that satisfy labor and land market
clearing conditions. The exogenous equilibrium is defined in page 
18 (2143) of the paper.  

If you are in the CLOSED-CITY case, then:

    a. H (total popolation) is an exogenous variable; and
    b. Ū (reservation/wider economy utility level) is an endogenous variable.

If you are in an OPEN-CITY case, then:

    a. H (total popolation) is an endogenous variable; and
    b. Ū (reservation/wider economy utility level) is an exogenous variable.

--- IGNORE ---

1.  This function solves for the equilibrium of the model by simultaneously
    iterating over guesses of w̃ⱼ, θᵢ and Qⱼ until convergence is achieved, 
    using a damping factor to ensure stability.
2.  In the ARSW tookit, the tolerance is implicitly of 2 digits. I will be 
    a little stricter and impose a 3-digit tolerance to leverage on Julia's
    higher computational efficiency. This is different from the rest of the
    functions which have a 6-digit tolerance as standard.
"""
function solve_equilibrium(params::ModelParameters, exo_fund::ExogenousFundamentals, pop_uti::Float64; endogenous_agglomeration::Bool=false, endo_params::Union{EndogenousModelParameters, Nothing}=nothing, prices_guess::Union{PricesGuess, Nothing} = nothing, tol_digits::Int=3, iter_max::Int=1000, damp_fact::Float64 = 0.4, open_city::Union{Bool,Nothing}=nothing)
    # checking if we are good to go
    if endogenous_agglomeration && isnothing(endo_params)
        ErrorException("The endogenous agglomeration equilibrium demands further the λ, δ, η, and ρ parameters.")
    elseif isnothing(open_city)
        ErrorException("Please specificy if you are in a open-city or closed-city model with the `open_city` keyword argument.")
    end

    # unpack parameters
    (; α, β, μ, κ, ε) = params; γ = gamma((ε-1)/ε);
    
    if endogenous_agglomeration
        # unpack remaining parameters
        (; λ, δ, η, ρ) = endo_params;
        # if in the endogenous agglomeration model, treat the first two entries of `exo_fund` as the exogenous part of the agglomeration forces
        aⱼ=copy(exo_fund.Ãⱼ); bᵢ=copy(exo_fund.B̃ᵢ); 
        φᵢ = copy(exo_fund.φᵢ); Kᵢ = copy(exo_fund.Kᵢ); τᵢⱼ = copy(exo_fund.τᵢⱼ);
    else
        # unpack exogenous fundamentals
        (; Ãⱼ, B̃ᵢ, φᵢ, Kᵢ, τᵢⱼ) = exo_fund;    
    end

    # separate last argument between the open- and closed-city cases
    Ū = open_city ? copy(pop_uti) : 1; Ūcity = copy(Ū); 
    H̃ = open_city ? (Ū/γ)^ε : copy(pop_uti); # initial guess w/ eq. 9

    # positional variables
    pos_employment = endogenous_agglomeration ? vec(aⱼ.>0) : vec(Ãⱼ.>0)
    pos_residence = endogenous_agglomeration ? vec(bᵢ.>0) : vec(B̃ᵢ.>0)
    idx_emp = findall(pos_employment) ; idx_res = findall(pos_residence);
    n_places = size(Kᵢ,1); n_workplaces = size(idx_emp,1); n_residence = size(idx_res,1);

    # initializing variables to be updated in the loop
    Qⱼ1 = zeros(n_places); w̃ⱼ1 = zeros(n_places);
    θᵢ1 = zeros(n_places); H̃ₘⱼ = zeros(n_places); 
    H̃ᵣᵢ = zeros(n_places); Ỹⱼ = zeros(n_places); 
    Ew̃ᵢ = zeros(n_places); Φᵢⱼ = zeros(n_residence, n_workplaces);
    Lₘⱼ = zeros(n_places); πᵢⱼ = zeros(n_places, n_places);

    # initial values for endogenous agglomeration forces; only relevant in the endogenous agglomeration model 
    if endogenous_agglomeration 
        # initial values for employement distribution (eq. 9)
        @. H̃ₘⱼ[pos_employment] =  H̃/n_workplaces;
        @. H̃ᵣᵢ[pos_residence] = H̃/n_residence;
        # initial values for productivities and amenities (eq. 20 + 21) 
        Ãⱼ = @. aⱼ * $sum(exp(-δ*τᵢⱼ)*(H̃ₘⱼ/Kᵢ)', dims=2)^λ ;
        B̃ᵢ = @. bᵢ * $sum(exp(-ρ*τᵢⱼ)*(H̃ᵣᵢ/Kᵢ)', dims=2)^η;
    end 
    
    # further positional arguments
    pure_emp = @. $vec((Ãⱼ>0) & (B̃ᵢ==0)); # I should use {H̃ₘⱼ,H̃ᵣᵢ} but it is identical to using {Ãⱼ,B̃ᵢ}
    pure_res = @. $vec((Ãⱼ==0) & (B̃ᵢ>0)); # I should use {H̃ₘⱼ,H̃ᵣᵢ} but it is identical to using {Ãⱼ,B̃ᵢ}
    shared_space = @. $vec((Ãⱼ>0) & (B̃ᵢ>0)); # I should use {H̃ₘⱼ,H̃ᵣᵢ} but it is identical to using {Ãⱼ,B̃ᵢ}

    # initial guess for the equilibrium prices of the model
    pg = isnothing(prices_guess) ? PricesGuess(n_places, pure_res, shared_space) : prices_guess
    Qⱼ0 = copy(pg.Qⱼ0); w̃ⱼ0 = copy(pg.w̃ⱼ0); θᵢ0 = copy(pg.θᵢ0);
    
    # completely specialized blocks never change since amenity or productivity are zero  
    @. θᵢ1[pure_emp] = 1;
    @. θᵢ1[pure_res] = 0;

    # defining loop variables
    iter = 0; tol = 10.0^(-tol_digits);
    err_Q = 1; err_w = 1; err_θ = 1; err_pop_uti = 1;

    # other variables 
    dᵢⱼ = @. exp(κ*τᵢⱼ[idx_res,idx_emp]); # by assumption
    Lᵢ = @. φᵢ * Kᵢ^(1-μ); # eq. 19

    # initiate the model loop
    if endogenous_agglomeration && open_city
        println(">>>> Solving the open-city endogenous agglomeration equilibrium <<<<") 
    elseif endogenous_agglomeration && !open_city  
        println(">>>> Solving the closed-city endogenous agglomeration equilibrium <<<<")
    elseif !endogenous_agglomeration && open_city
        println(">>>> Solving the open-city exogenous agglomeration equilibrium <<<<")
    else 
        println(">>>> Solving the closed-city exogenous agglomeration equilibrium <<<<")
    end
    
    while ((err_Q ≥ tol) || (err_w ≥ tol) || (err_θ ≥ tol) || (err_pop_uti ≥ tol)) && (iter ≤ iter_max)
        # updating endogenous variables by solving the model equations

        # --- πᵢⱼ through eq. 4 ---
        @. Φᵢⱼ = (B̃ᵢ[pos_residence] * w̃ⱼ0[pos_employment]')^ε * (dᵢⱼ*Qⱼ0[pos_residence]^(1-β))^(-ε);
        Φ =  sum(Φᵢⱼ); # Φ need not match H̃; only a specific rescaling of B̃ᵢ makes it coincide, so it diverges in counterfactual/comparative statics exercises.
        @. πᵢⱼ[idx_res,idx_emp] = Φᵢⱼ / Φ;

        # --- H̃ₘⱼ and H̃ᵣᵢ through eq. 5 ---
        @. H̃ₘⱼ = $sum(πᵢⱼ, dims=1)' * H̃ ;
        @. H̃ᵣᵢ = $sum(πᵢⱼ, dims=2) * H̃ ;
        
        # --- spatial equilibrium condition through eq. 9 ---
        Ūcity = γ * Φ ^ (1/ε);

        # --- If in the open-city case, adjust population ---
        open_city && (H̃1 = (Ūcity/Ū)^ε * H̃); # increase population if within-city utility exceeds that of the wider economy; from eq. 9 one can also infer that employment scales in utility at elasticity ε. 

        # --- If in the endogenous agglomeration case, update Ãⱼ and B̃ᵢ through eq. 20 + 21 ---
        if endogenous_agglomeration
            @. Ãⱼ = aⱼ * $sum(exp(-δ*τᵢⱼ)*(H̃ₘⱼ/Kᵢ)', dims=2)^λ ; 
            @. B̃ᵢ = bᵢ * $sum(exp(-ρ*τᵢⱼ)*(H̃ᵣᵢ/Kᵢ)', dims=2)^η;
        end

        # --- w̃ⱼ through eq. 10 + eq. 11 ---
        @. Ỹⱼ = Ãⱼ * H̃ₘⱼ^α * (θᵢ0 * Lᵢ)^(1-α); 
        @. w̃ⱼ1[pos_employment] = α * Ỹⱼ[pos_employment] / H̃ₘⱼ[pos_employment];

        # --- Qᵢ thorugh eq. S.20 + eq. 17 + eq. 18 + eq. 14 ---
        @. Ew̃ᵢ[pos_residence] = $sum(πᵢⱼ[idx_res,idx_emp] / $sum(πᵢⱼ[idx_res,idx_emp], dims=2) * w̃ⱼ0[pos_employment]' , dims=2); 
        @. Qⱼ1[pure_res] = ((1-β) * Ew̃ᵢ[pure_res] * H̃ᵣᵢ[pure_res]) / ((1-θᵢ0[pure_res]) * Lᵢ[pure_res]);
        @. Qⱼ1[pure_emp] = ((1-α) * Ỹⱼ[pure_emp]) / (θᵢ0[pure_emp] * Lᵢ[pure_emp]); # you could maybe use S.24 here?!
        @. Qⱼ1[shared_space] = (((1-α) * Ỹⱼ[shared_space]) + ((1-β) * Ew̃ᵢ[shared_space] * H̃ᵣᵢ[shared_space])) / Lᵢ[shared_space]; # akin to: θ * Qⱼ1[pure_emp] + (1-θ) * Qⱼ1[pure_res]

        # --- θᵢ through eq. 10 (FOC wrt Lₘ, i.e., S.23) + S.53 --- PS.: my understanding is that since Lᵢ is exogenous, S.53 guarantees a market-clearing equilibrium.
        @. Lₘⱼ[pos_employment] = (1-α) * Ỹⱼ[pos_employment] / Qⱼ0[pos_employment]; # why could you not use S.49 here? Maybe because I have H̃ⱼ instead of Hⱼ... 
        @. θᵢ1[shared_space] = Lₘⱼ[shared_space] / (Lᵢ[shared_space]);

        # update error metrics here
        iter += 1; 
        err_Q = @. $maximum(abs(Qⱼ1 - Qⱼ0)); 
        err_w = @. $maximum(abs(w̃ⱼ1 - w̃ⱼ0)); 
        err_θ = @. $maximum(abs(θᵢ1 - θᵢ0)); 
        open_city ? (err_pop_uti = abs(Ūcity-Ū)) : (err_pop_uti = 0); # not checking for convergency, but for consistency w/ the data

        # revise guesses
        @. Qⱼ0 = (1-damp_fact) * Qⱼ0 + damp_fact * Qⱼ1 ;
        @. w̃ⱼ0 = (1-damp_fact) * w̃ⱼ0 + damp_fact * w̃ⱼ1 ;
        @. θᵢ0 = (1-damp_fact) * θᵢ0 + damp_fact * θᵢ1 ;
        open_city && (H̃ = (1-(0.2*damp_fact)) * H̃ + (0.2*damp_fact) * H̃1) ; # only in the open city case. Damping the `damp_fact` even further because there is a strong positive feedback loop here, which is amplified due to ε being relatively large 

        # Print convergence rate
        println([iter, trunc(err_Q / tol, digits=0), trunc(err_w / tol, digits=0), trunc(err_θ / tol, digits=0), trunc(err_pop_uti / tol, digits=0)]);
    end
    
    # Print status
    (iter < iter_max) ? println(">>>> Equilibrium achieved! <<<<") : println(">>>> Failed to find an equilibrium <<<<")

    # Return the equilibrium endogenous variables 
    return open_city ? (Qⱼ0, w̃ⱼ0, θᵢ0, πᵢⱼ, H̃) : (Qⱼ0, w̃ⱼ0, θᵢ0, πᵢⱼ, Ūcity) 
end

"""
The function detangles the overall agglomeration forces {Ãⱼ,B̃ᵢ} into
its endogenous and exogenous components {Υⱼ,Ωᵢ} and {aⱼ, bᵢ}, respectivelly.
It assumes you provide the properly estimated {λ, δ, η, ρ} set of parameters.
Importantly, it only reports the exogenous component of the agglomeration
forces, since it is the only part that matters for estimating the GE model.   
"""
function detangle_agglomeration(params::EndogenousModelParameters, exo_fund::ExogenousFundamentals, Hₘⱼ::Vector{Float64}, Hᵣᵢ::Vector{Float64})
    # unpacking parameters
    (; λ, δ, η, ρ) = params;
    (; Ãⱼ, B̃ᵢ, φᵢ, Kᵢ, τᵢⱼ) = exo_fund;

    # Detangling productivity with equation (20)
    Υⱼ = @. $sum(exp(-δ*τᵢⱼ)*(Hₘⱼ/Kᵢ)', dims=2);
    aⱼ = @. Ãⱼ / (Υⱼ^λ);

    # Detangling amenities with equation (21)
    Ωᵢ = @. $sum(exp(-ρ*τᵢⱼ)*(Hᵣᵢ/Kᵢ)', dims=2);
    bᵢ = @. B̃ᵢ / (Ωᵢ^η);

    return vec(aⱼ), vec(bᵢ)
end

end
