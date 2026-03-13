function cal_model_seq(Qⱼ,Hₘⱼ,Hᵣᵢ,τᵢⱼ,Kᵢ; tol_digits=6)
    "
    This function assumes that you have predefined the parameters
    ε, κ, α, β, and μ. It then computes the structural fundamentals 
    of the model given the observed ('real world') equilibrium
    variables:
        1. Qⱼ = rent prices; 
        2. Hₘⱼ = workplace employment (population);
        3. Hᵣᵢ = residential employment (population);
    And the exogenous fundamentals
        1. τᵢⱼ = bilateral travel time matrix s.t. rows (i) denote 
            residences and columns (j) denote workplaces; and
        2. Kᵢ = geographical area (unkown unit).
    The output of this function is the set of structural fundamentals and
    endogenous variales of of the model (Ãⱼ, B̃ᵢ, w̃ⱼ, πᵢⱼ, Tw̃ᵢ, ϕᵢ, Lᵢᴰ, 
    θᵢ, H̃ₘⱼ, H̃ᵣᵢ, CMA) that are consistent with the observed equilibrium.
    --- IGNORE ---
    This function solves for the equilibrium of the model by iterating 
    over wages to assess for the equilibrium productivity. All else is
    sequentially (and algebraically) derived from these results, being 
    (re)scaled to match the data.
    "
    # Identifying places with firms and residents
    pos_employment = vec(Hₘⱼ.>0); pos_residence = vec(Hᵣᵢ.>0); 
    idx_emp = findall(pos_employment); idx_res = findall(pos_residence);
    n_places = size(Qⱼ,1);

    # **************************************************
    # *** w̃ⱼ (adjusted wages) and Ãⱼ (productivity) ****
    # **************************************************
    
    # computing transformed wages (ωⱼ) array
    ωⱼ, Ĥₘⱼ = get_ω(Hₘⱼ,Hᵣᵢ,τᵢⱼ,Qⱼ, tol_digits=tol_digits, ε=ε); 
    w̃ⱼ = @. ωⱼ ^ (1 / ε);  # recover adjusted wages by remembering that w̃ⱼ = ω^(1/ε) = wⱼEⱼ^(1/ε)
    @. w̃ⱼ[w̃ⱼ .> 0] = w̃ⱼ[w̃ⱼ .> 0] / $geomean(w̃ⱼ[w̃ⱼ .> 0]); # normalizing adjusted wages
    
    # Compute adjusted productivity (from equation 12) up to scale (due to wages)
    Ãⱼ = @. ((Qⱼ / (1 - α)) ^ (1 - α)) * ((w̃ⱼ / α) ^ α);
    
    # *******************************************************
    # *** B̃ᵢ (adjusted amenities) and CMA (Market Access) ***
    # *******************************************************

    # Commuting market access (CMA) from eq. (29)
    dᵢⱼ= @. exp(κ * τᵢⱼ[idx_res,idx_emp]); # iceberg commuting cost, by assumption
    CMA = zeros(n_places);
    @. CMA[pos_residence] = $sum(ωⱼ[pos_employment]' / (dᵢⱼ^ε), dims=2); 

    # Amenities from equation (28) or (S.47)
    B̃ᵢ = zeros(n_places); 
    Hᵣᵢᵃ = @. Hᵣᵢ[pos_residence] / $geomean(Hᵣᵢ[pos_residence]);
    Qⱼᵃ = @. Qⱼ[pos_residence] / $geomean(Qⱼ[pos_residence]);
    CMAᵃ = @. CMA[pos_residence] / $geomean(CMA[pos_residence]); 
    @. B̃ᵢ[pos_residence] = (Hᵣᵢᵃ)^(1/ε) * (Qⱼᵃ)^(1-β) * (CMAᵃ)^(-1/ε) ;
    
    # *******************************************************************
    # *** Rescaling Ãⱼ, B̃ᵢ, and computing  πᵢⱼ (commuting flow prob.) ***
    # *******************************************************************
        
    # Normalize productivity to geomean 1
    @. Ãⱼ[pos_employment] = Ãⱼ[pos_employment] / $geomean(Ãⱼ[pos_employment]);
 
    # Change wages and CMA to be consistent with the normalization on productivity (eq. 12)
    @. w̃ⱼ[pos_employment] = (Ãⱼ[pos_employment]^(1/α))*α*((1-α)/Qⱼ[pos_employment])^((1-α)/α);
    @. CMA[pos_residence] = $sum((w̃ⱼ[pos_employment]'/dᵢⱼ)^ε, dims=2);

    # Compute bilateral commuting probabilities (eq. 4)
    πᵢⱼ = zeros(n_places,n_places);
    Φᵢⱼ = @. (B̃ᵢ[pos_residence]*w̃ⱼ[pos_employment]')^ε * (dᵢⱼ*Qⱼ[pos_residence]^(1-β))^(-ε); # total population in the model
    @. πᵢⱼ[idx_res,idx_emp] = Φᵢⱼ / $sum(Φᵢⱼ);

    # Normalizing amenities to match data population
    @. B̃ᵢ[pos_residence] = B̃ᵢ[pos_residence] * ($sum(Hₘⱼ) / $sum(Φᵢⱼ))^(1/ε);
    "
    The authors measure utility in a unit measure s.t. (Ū/γ)ᵋ/H = 1, where γ = Γ(ε−1/ε) and Γ(·) is the Gamma function (See supplement p. 17).
    Thus, it is implied that Φ = H, as demonstrated in p. 18 of the supplement. Hence, if the population in the data (H) is greater than the 
    population in the model (Φ), we increase the amenities to make the city more attractive and attract more residents.
    "

    # ******************************************************
    # *** Tw̃ᵢ (total expected worker residential income) *** 
    # ******************************************************

    # Residential and Workplace probabilities (equation 5)
    πᵣᵢ = sum(πᵢⱼ, dims=2);
    πₘⱼ = sum(πᵢⱼ, dims=1)';

    # Predicted residence and workplace employment
    H̃ₘⱼ = @. πₘⱼ * $sum(Hₘⱼ);
    H̃ᵣᵢ = @. πᵣᵢ * $sum(Hᵣᵢ);

    # Compute expected residential work income (eq. S20)
    Ew̃ᵢ = zeros(n_places);
    @. Ew̃ᵢ[pos_residence] = $sum(πᵢⱼ[idx_res,idx_emp] / πᵣᵢ[pos_residence] * w̃ⱼ[pos_employment]' , dims=2);

    # Compute total expected residential worker income
    Tw̃ᵢ = @. Ew̃ᵢ * H̃ᵣᵢ;

    # ******************************
    # *** Density of development *** 
    # ******************************

    # Compute commercial/workplace floorspace demand (equation 18/S30)
    Lᵢᴹ = @. ((1-α)* Ãⱼ / Qⱼ)^(1/α) * Hₘⱼ;

    # Compute residential floorspace demand (equation 17/S29)
    Lᵢᴿ = @. (1-β) * Tw̃ᵢ / Qⱼ ;

    # Total floor space demand (by definition)
    Lᵢᴰ = @. Lᵢᴹ + Lᵢᴿ ;

    # Share of commercial floor space (by definition)
    θᵢ = @. Lᵢᴹ / Lᵢᴰ;

    # Density of development (equation 19/S31)
    ϕᵢ = @. Lᵢᴰ/(Kᵢ^(1-μ));
    
    return Ãⱼ, B̃ᵢ, w̃ⱼ, πᵢⱼ, Tw̃ᵢ, ϕᵢ, Lᵢᴰ, θᵢ, H̃ₘⱼ, H̃ᵣᵢ, CMA
end

function cal_model_sim(Qⱼ,Hₘⱼ,Hᵣᵢ,τᵢⱼ,Kᵢ; tol_digits=6, iter_max=1000)
    "
    This function assumes that you have predefined the parameters
    ε, κ, α, β, and μ. It then computes the structural fundamentals 
    of the model given the observed ('real world') equilibrium
    variables:
        1. Qⱼ = rent prices; 
        2. Hₘⱼ = workplace employment (population);
        3. Hᵣᵢ = residential employment (population);
    and structural fundamentals:
        1. τᵢⱼ = bilateral travel time matrix s.t. rows (i) denote 
            residences and columns (j) denote workplaces; and
        2. Kᵢ = geographical area (unkown unit).
    The output of this function is the set of structural fundamentals and
    endogenous variables of the model (Ãⱼ, B̃ᵢ, w̃ⱼ, πᵢⱼ, Tw̃ᵢ, ϕᵢ, Lᵢᴰ, θᵢ, 
    H̃ₘⱼ, H̃ᵣᵢ, CMA) that are consistent with the observed equilibrium.
    --- IGNORE ---
    This function solves for the equilibrium of the model by simultaneously
    iterating over guesses of Ãⱼ and B̃ᵢ up until convergence is achieved.
    "
    # Identifying places with firms and residents
    pos_employment = vec(Hₘⱼ.>0); pos_residence = vec(Hᵣᵢ.>0); 
    idx_emp = findall(pos_employment); idx_res = findall(pos_residence);
    n_places = size(Qⱼ,1);

    # Defining initial guesses
    Ãⱼ0 = zeros(n_places); Ãⱼ0[pos_employment] .= 1;
    B̃ᵢ0 = zeros(n_places); B̃ᵢ0[pos_residence] .= 1;

    # initiating variables to be updated in the loop
    w̃ⱼ = zeros(n_places); πᵢⱼ = zeros(n_places,n_places);
    Tw̃ᵢ = zeros(n_places); Lᵢᴿ = zeros(n_places);
    Lᵢᴹ = zeros(n_places); Ãⱼ1 = zeros(n_places); 
    B̃ᵢ1 = zeros(n_places); H̃ₘⱼ = zeros(n_places); 
    H̃ᵣᵢ = zeros(n_places); Φᵢⱼ = zeros(size(idx_res,1),size(idx_emp,1));
    CMA = zeros(n_places);
    
    # Setting up convergence criteria and additional variables
    iter = 0; err_Ãⱼ = 10000; err_B̃ᵢ = 10000; tol = 10.0^(-tol_digits);
    dᵢⱼ = @. exp(κ*τᵢⱼ[idx_res,idx_emp]); # iceberg commuting cost, by assumption
    
    # initiate the model loop
    println(">>>> Calibrating Ã and B̃ <<<<")
    while  ((err_Ãⱼ >= tol) || (err_B̃ᵢ >= tol)) && (iter <= iter_max)
        
        # Guess wages using the first-order condition (eq. 12)
        @. w̃ⱼ[pos_employment] = (((1-α)/Qⱼ[pos_employment])^((1-α)/α))*α*(Ãⱼ0[pos_employment]^(1/α));
        
        # Compute bilateral commuting probabilities (eq. 4)
        @. Φᵢⱼ = (B̃ᵢ0[pos_residence]*w̃ⱼ[pos_employment]')^ε * (dᵢⱼ*Qⱼ[pos_residence]^(1-β))^(-ε); # total population in the model
        @. πᵢⱼ[idx_res,idx_emp] = Φᵢⱼ / $sum(Φᵢⱼ); # unconditional commuting probabilities

        # Compute predicted residence and workplace employment (eq. 5)
        @. H̃ₘⱼ = $sum(πᵢⱼ, dims=1)' * $sum(Hₘⱼ);
        @. H̃ᵣᵢ = $sum(πᵢⱼ, dims=2) * $sum(Hₘⱼ);

        # Updating guesses
        @. Ãⱼ1[pos_employment] = (Hₘⱼ[pos_employment]/H̃ₘⱼ[pos_employment])^(1/ε) * Ãⱼ0[pos_employment]; # slightly increase productivity if predicted employment is lower than data
        @. B̃ᵢ1[pos_residence] = (Hᵣᵢ[pos_residence]/H̃ᵣᵢ[pos_residence])^(1/ε) * B̃ᵢ0[pos_residence]; # slightly increase amenities if predicted population is lower than data
        
        # Check if updated values are valid (i.e. non-nan)
        if (sum(isnan.(Ãⱼ1)) > 0) || (sum(isnan.(B̃ᵢ1)) > 0)
            # set to random values around 1
            @. Ãⱼ1[pos_employment] = 0.95 + (1.05-0.95) * rand($length(Ãⱼ1[pos_employment]));
            @. B̃ᵢ1[pos_residence] = 0.95 + (1.05-0.95) * rand($length(B̃ᵢ1[pos_residence]));
        end
        
        # Damping the updates to improve stability (I will follow ARSW toolkit and use a 0.5 damping factor, even if 0.75/0.25 should be safer)
        @. Ãⱼ0 = 0.5 * Ãⱼ0 + 0.5 * Ãⱼ1 ;
        @. B̃ᵢ0 = 0.5 * B̃ᵢ0 + 0.5 * B̃ᵢ1 ;

        # Normalizing productivity to geomean 1
        @. Ãⱼ0[pos_employment] = Ãⱼ0[pos_employment] / $geomean(Ãⱼ0[pos_employment]);

        # Normalizing amenities to match data population
        @. B̃ᵢ0[pos_residence] = B̃ᵢ0[pos_residence] * ($sum(Hₘⱼ)/$sum(Φᵢⱼ))^(1/ε);

        # Update iteration variables
        iter += 1; 
        err_Ãⱼ = @. $round($maximum(abs(Ãⱼ1 - Ãⱼ0)),digits=tol_digits); 
        err_B̃ᵢ = @. $round($maximum(abs(B̃ᵢ1 - B̃ᵢ0)),digits=tol_digits);

        # Print convergence rate
        println([iter, trunc(err_Ãⱼ / tol, digits=0), trunc(err_B̃ᵢ / tol, digits=0)])
    end
    if iter==iter_max
        error("Convergence not achieved for adjusted wages Ã and B̃")
    end
    println(">>>> Ã and B̃ Converged <<<<")

    # Compute total expected residential worker income (eq. S20)
    @. Tw̃ᵢ[pos_residence] = $sum(πᵢⱼ[idx_res,idx_emp] / $sum(πᵢⱼ, dims=2)[pos_residence] * w̃ⱼ[pos_employment]' , dims=2) * H̃ᵣᵢ[pos_residence];

    # Compute CMA (equation 29)
    @. CMA[pos_residence] = $sum((w̃ⱼ[pos_employment]' / dᵢⱼ) ^ ε , dims=2);

    # Compute residential/commertial floorspace (equations S29 and S30)
    @. Lᵢᴿ[pos_residence] = (1-β) * Tw̃ᵢ[pos_residence] / Qⱼ[pos_residence];
    @. Lᵢᴹ[pos_employment] = ((1-α) * Ãⱼ0[pos_employment] / Qⱼ[pos_employment])^(1/α) * H̃ₘⱼ[pos_employment];
    Lᵢᴰ = @. Lᵢᴿ + Lᵢᴹ; 
    
    # Compute density of development ϕᵢ (equation S.31)
    ϕᵢ = @. Lᵢᴰ / (Kᵢ ^ (1-μ));

    # Compute commercial floor space share θᵢ (definition) 
    θᵢ = @. Lᵢᴹ / Lᵢᴰ;

    return Ãⱼ0, B̃ᵢ0, w̃ⱼ, πᵢⱼ, Tw̃ᵢ, ϕᵢ, Lᵢᴰ, θᵢ, H̃ₘⱼ, H̃ᵣᵢ, CMA
end 