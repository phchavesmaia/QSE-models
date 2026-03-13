module ModelSolver

using SpecialFunctions
export ModelParameters, ExogenousFundamentals, PricesGuess, solve_equilibrium
struct ModelParameters
    α::Real # Share Of consumption
    β::Real # Share Of Labor Costs
    κ::Real # Size Of Commuting Costs
    ε::Real # Fréchet Shape Parameter
    μ::Real # Share Of Capital In Construction Costs
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

function PricesGuess(n_places::Int, pure_res, shared_space)
    Q = ones(n_places)
    w = ones(n_places)
    θ = ones(n_places)
    θ[pure_res] .= 0
    θ[shared_space] .= 0.5
    return PricesGuess(Q, w, θ)
end

function solve_equilibrium(params::ModelParameters, exo_fund::ExogenousFundamentals, pop_uti::Real; prices_guess::Union{PricesGuess, Nothing} = nothing, tol_digits=3, iter_max=1000, damp_fact = 0.5, closed_city=true)
    "
        This function assumes that you have predefined the parameters
        {α, β, κ, ε, and μ}. It then solves for the general equilibrium 
        of the model given the structural exogenous fundamentals:
            1. Ãⱼ = workplace productivity;
            2. B̃ᵢ = residential amenities;
            3. φᵢ = density of development;
            4. Kᵢ = geographical area; and
            5. τᵢⱼ = bilateral travel time matrix.
        The output of this function is the set of endogenous equilibrium
        variables {w̃ⱼ, θᵢ, Qⱼ, πᵢⱼ} that satisfy labor and land market
        clearing conditions. The equilibrium is defined in page 18 (2143)
        of the paper. 

        If we are in the CLOSED-CITY case, then
            a. H (total popolation) is an exogenous variable; and
            b. Ū (reservation/wider economy utility level) is an endogenous variable.
        If we are in an OPEN-CITY case, then
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
    "
    # unpack parameters
    (; α, β, κ, ε, μ) = params; γ = gamma((ε-1)/ε);
    (; Ãⱼ, B̃ᵢ, φᵢ, Kᵢ, τᵢⱼ) = exo_fund;

    # separate last argument between the open- and closed-city cases
    H̃ = closed_city ? pop_uti : 1;
    Ū = closed_city ? 1 : pop_uti;   

    # positional variables
    pos_employment = vec(Ãⱼ.>0) ; pos_residence = vec(B̃ᵢ.>0);
    idx_emp = findall(pos_employment) ; idx_res = findall(pos_residence);
    n_places = size(Kᵢ,1); n_workplaces = size(idx_emp,1); n_residence = size(idx_res,1);
    pure_emp = @. (Ãⱼ>0) & (B̃ᵢ==0); # I should use {H̃ₘⱼ,H̃ᵣᵢ} but it is identical to using {Ãⱼ,B̃ᵢ}
    pure_res = @. (Ãⱼ==0) & (B̃ᵢ>0); # I should use {H̃ₘⱼ,H̃ᵣᵢ} but it is identical to using {Ãⱼ,B̃ᵢ}
    shared_space = @. (Ãⱼ>0) & (B̃ᵢ>0); # I should use {H̃ₘⱼ,H̃ᵣᵢ} but it is identical to using {Ãⱼ,B̃ᵢ}

    # initial guess for the equilibrium prices of the model
    (; Qⱼ0, w̃ⱼ0, θᵢ0) = isnothing(prices_guess) ? PricesGuess(n_places, pure_res, shared_space) : prices_guess

    # initializing variables to be updated in the loop
    Qⱼ1 = zeros(n_places); w̃ⱼ1 = zeros(n_places);
    θᵢ1 = zeros(n_places); H̃ₘⱼ = zeros(n_places); 
    H̃ᵣᵢ = zeros(n_places); Ỹⱼ = zeros(n_places); 
    Ew̃ᵢ = zeros(n_places); Φᵢⱼ = zeros(n_residence, n_workplaces);
    Lₘⱼ = zeros(n_places); πᵢⱼ = zeros(n_places, n_places);

    # completely specialized blocks never change since amenity or productivity are zero  
    @. θᵢ1[pure_emp] = 1;
    @. θᵢ1[pure_res] = 0;
    
    # defining loop variables
    iter = 0; tol = 10.0^(-tol_digits);
    err_Q = 10000; err_w = 10000; err_θ = 10000;
    
    # other variables 
    dᵢⱼ = @. exp(κ*τᵢⱼ[idx_res,idx_emp]); # by assumption
    Lᵢ = @. φᵢ * Kᵢ^(1-μ); # eq. 19

    # initiate the model loop
    if closed_city
        println(">>>> Solving the closed-city equilibrium <<<<")
    else
        println(">>>> Solving the open-city equilibrium <<<<")
    end

    while ((err_Q >= tol) || (err_w >= tol) || (err_θ >= tol)) && (iter <= iter_max)
        # updating endogenous variables by solving the model equations

        # --- πᵢⱼ through eq. 4 ---
        @. Φᵢⱼ = (B̃ᵢ[pos_residence] * w̃ⱼ0[pos_employment]')^ε * (dᵢⱼ*Qⱼ0[pos_residence]^(1-β))^(-ε);
        Φ =  sum(Φᵢⱼ); 
        @. πᵢⱼ[idx_res,idx_emp] = Φᵢⱼ / Φ;
        
        # --- spatial equilibrium condition through eq. 9 ---
        if closed_city
            # Ū is endogenous in the closed-city equilibrium
            Ū = γ * Φ ^ (1/ε);
        else
            # H̃ is endogenous in the open-city equilibrium
            Ūcity = γ * Φ ^ (1/ε);
            H̃ = (Ūcity/Ū)^ε * H̃; # increase population if within-city utility exceeds that of the wider economy; from eq. 9 one can also infer that employment scales in utility at elasticity ε. 
        end

        # --- H̃ₘⱼ and H̃ᵣᵢ through eq. 5 ---
        @. H̃ₘⱼ = $sum(πᵢⱼ, dims=1)' * H̃ ;
        @. H̃ᵣᵢ = $sum(πᵢⱼ, dims=2) * H̃ ;
        
        # --- w̃ⱼ through eq. 10 + eq. 11 ---
        @. Ỹⱼ = Ãⱼ * H̃ₘⱼ^α * (θᵢ0 * Lᵢ)^(1-α); 
        @. w̃ⱼ1[pos_employment] = α * Ỹⱼ[pos_employment] / H̃ₘⱼ[pos_employment];

        # --- Qᵢ thorugh eq. S.20 + eq. 17 + eq. 18 + eq. 14 ---
        @. Ew̃ᵢ[pos_residence] = $sum(πᵢⱼ[idx_res,idx_emp] / $sum(πᵢⱼ[idx_res,idx_emp], dims=2) * w̃ⱼ0[pos_employment]' , dims=2); 
        @. Qⱼ1[pure_res] = ((1-β) * Ew̃ᵢ[pure_res] * H̃ᵣᵢ[pure_res]) / ((1-θᵢ0[pure_res]) * Lᵢ[pure_res]);
        @. Qⱼ1[pure_emp] = ((1-α) * Ỹⱼ[pure_emp]) / (θᵢ0[pure_emp] * Lᵢ[pure_emp]); # you could maybe use S.24 here?!
        @. Qⱼ1[shared_space] = (((1-α) * Ỹⱼ[shared_space]) + ((1-β) * Ew̃ᵢ[shared_space] * H̃ᵣᵢ[shared_space])) / Lᵢ[shared_space]; # akin to: θ * Qⱼ1[pure_emp] + (1-θ) * Qⱼ1[pure_res]

        # --- θᵢ through eq. 10 (CPO wrt Lₘ, i.e., S.23) + S.53 --- PS.: my understanding is that since Lᵢ is exogenous, S.53 guarantees a market-clearing equilibrium.
        @. Lₘⱼ[pos_employment] = (1-α) * Ỹⱼ[pos_employment] / Qⱼ0[pos_employment]; # you could maybe use S.49 here?!
        @. θᵢ1[shared_space] = Lₘⱼ[shared_space] / (Lᵢ[shared_space]);

        # update error metrics here
        iter += 1; 
        err_Q = @. $round($maximum(abs(Qⱼ1 - Qⱼ0)),digits=tol_digits); 
        err_w = @. $round($maximum(abs(w̃ⱼ1 - w̃ⱼ0)),digits=tol_digits); 
        err_θ = @. $round($maximum(abs(θᵢ1 - θᵢ0)),digits=tol_digits); 

        # revise guesses (safer damping; otherwise it 'bounces' a lot)
        @. Qⱼ0 = (1-damp_fact) * Qⱼ0 + damp_fact * Qⱼ1 ;
        @. w̃ⱼ0 = (1-damp_fact) * w̃ⱼ0 + damp_fact * w̃ⱼ1 ;
        @. θᵢ0 = (1-damp_fact) * θᵢ0 + damp_fact * θᵢ1 ;

        # Print convergence rate
        println([iter, trunc(err_Q / tol, digits=0), trunc(err_w / tol, digits=0), trunc(err_θ / tol, digits=0)])
    end
    
    # Print status
    if iter < iter_max
        println(">>>> Equilibrium achieved! <<<<")
    else
        println(">>>> Failed to find an equilibrium <<<<")
    end

    # Return the equilibrium endogenous variables 
    if closed_city
        return Qⱼ0, w̃ⱼ0, θᵢ0, πᵢⱼ, Ū
    else
        return Qⱼ0, w̃ⱼ0, θᵢ0, πᵢⱼ, H̃
    end
end

end