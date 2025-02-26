function cal_model(Qⱼ,Hₘⱼ,Hᵣᵢ,τᵢⱼ,Kᵢ; tol_digits=6)

    # Identifying places with firms and residents
    pos_employment = vec(Hₘⱼ.>0); pos_residence = vec(Hᵣᵢ.>0) 

    # **************************************************
    # *** w̃ⱼ (transformed wages) and Ãⱼ productivity ***
    # **************************************************
    
    # computing adjusted wages (ωⱼ) array
    ωⱼ, Ĥₘⱼ = get_ω(Hₘⱼ,Hᵣᵢ,τᵢⱼ,Qⱼ, tol_digits=tol_digits); 
    w̃ⱼ = ωⱼ .^ (1 / ε);  # remember that w̃ⱼ = ω^(1/ε) = wⱼEⱼ^(1/ε)
    w̃ⱼ[w̃ⱼ .> 0] = w̃ⱼ[w̃ⱼ .> 0] ./ geomean(w̃ⱼ[w̃ⱼ .> 0]); # normalizing transformed wages

    # 'validating' estimates
    df = DataFrame();
    df[!,"model"] = vec(Ĥₘⱼ);
    df[!,"data"] = vec(Hₘⱼ);
    model = reg(df,@formula(data ~ model));
    slope = round(coef(model)[2],digits=tol_digits);
    println("The slope of the model-real workplace population data is: $slope") # should be 1!
    
    # Compute adjusted productivity (from equation 12) up to scale (due to wages)
    Ãⱼ = ((Qⱼ ./ (1 - α)) .^ (1 - α)) .* ((w̃ⱼ ./ α) .^ α);
    
    # **********************************************
    # *** B̃ᵢ (amenities) and Market Access (CMA) ***
    # **********************************************

    # Commuting market access (CMA) from eq. (29)
    CMA = sum(ωⱼ'./exp.(ν.*τᵢⱼ'), dims=2); CMAₐ = CMA[pos_residence]; CMAₐ = CMAₐ./geomean(CMAₐ); 

    # Amenities from equation (28) or (S.47)
    B̃ᵢ = zeros(size(Qⱼ,1),1); 
    B̃ᵢ[pos_residence] = (Hᵣᵢ[pos_residence]./geomean(Hᵣᵢ[pos_residence])).^(1/ε) .*   (Qⱼ[pos_residence]./geomean(Qⱼ[pos_residence])).^(1-β) .* (CMAₐ).^(-1/ε) ;
    
    # *******************************************************************
    # *** Rescaling Ãⱼ, B̃ᵢ, and computing  πᵢⱼ (commuting flow prob.) ***
    # *******************************************************************

    # Normalize productivity to geomean 1
    Ãⱼ[pos_employment] = Ãⱼ[pos_employment]./geomean(Ãⱼ[pos_employment])
 
    # Change wages to be consistent the normalization on productivity (eq. 12)
    w̃ⱼ[pos_employment] = (Ãⱼ[pos_employment].^(1/α)).*α.*((1-α)./Qⱼ[pos_employment]).^((1-α)/α)

    # Compute bilateral commuting probabilities (eq. 4)
    πᵢⱼ = zeros(size(Hᵣᵢ,1),size(Hₘⱼ,1)); dᵢⱼ= exp.(κ.*τᵢⱼ[findall(pos_employment),findall(pos_residence)]')
    ϕᵢⱼ = (B̃ᵢ[pos_residence].*w̃ⱼ[pos_employment]').^ε .* (dᵢⱼ.*Qⱼ[pos_residence].^(1-β)).^(-ε);	
    πᵢⱼ[findall(pos_residence),findall(pos_employment)] = ϕᵢⱼ ./ sum(ϕᵢⱼ);

    # Normalizing amenities to match data population
    B̃ᵢ[pos_residence] = B̃ᵢ[pos_residence] .* (sum(Hₘⱼ)./sum(ϕᵢⱼ)).^(1 ./ ε)
    "
    The authors measure utility in a unit measure s.t. (Ū/γ)ᵋ/H = 1, where γ = Γ(ε−1/ε) and Γ(·) is the Gamma function (See supplement p. 17).
    Thus, it is implied that ϕ = H, as demonstrated in p. 18 of the supplement. Thus, if the population in the data (H) is greater than the 
    population in the model (ϕ), we increase the amenities to make the city more attractive and attract more residents.
    "

    # ******************************************************
    # *** Tw̃ᵢ (total expected worker residential income) *** 
    # ******************************************************

    # Residential and Workplace probabilities (equation 5)
    πᵣᵢ = sum(πᵢⱼ, dims=2);
    πₘⱼ = sum(πᵢⱼ, dims=1)';

    # Predicted residence and workplace employment
    H̃ₘⱼ = πₘⱼ .* sum(Hₘⱼ);
    H̃ᵣᵢ = πᵣᵢ .* sum(Hₘⱼ);

    # Compute expected residential work income (eq. S20)
    Ew̃ᵢ = zeros(size(Hᵣᵢ,1),1);
    Ew̃ᵢ[pos_residence] = sum(πᵢⱼ[findall(pos_residence),findall(pos_employment)] ./ πᵣᵢ[pos_residence] .* w̃ⱼ[pos_employment]' , dims=2);

    # Compute total expected residential worker income
    Tw̃ᵢ = Ew̃ᵢ .* H̃ᵣᵢ;

    # ******************************
    # *** Density of development *** 
    # ******************************

    # Compute commercial/workplace floorspace demand (equation 18/S30)
    Lᵢᴹ = ((1-α).* Ãⱼ ./ Qⱼ).^(1/α) .* Hₘⱼ;

    # Compute residential floorspace demand (equation 17/S29)
    Lᵢᴿ = (1-β) .* Tw̃ᵢ ./ Qⱼ ;

    # Total floor space demand (by definition)
    Lᵢᴰ = Lᵢᴹ + Lᵢᴿ ;

    # Share of commercial floor space (by definition)
    θᵢ = Lᵢᴹ ./ Lᵢᴰ;

    # Density of development (equation 19/S31)
    ϕᵢ = Lᵢᴰ./(Kᵢ.^(1-μ));

    return Ãⱼ, B̃ᵢ, w̃ⱼ, πᵢⱼ, Tw̃ᵢ, ϕᵢ, Lᵢᴰ, θᵢ, H̃ₘⱼ, H̃ᵣᵢ, CMA
end
