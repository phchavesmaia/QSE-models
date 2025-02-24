function cal_exog(Qⱼ,Hₘⱼ,Hᵣᵢ,τᵢⱼ; tol_digits=6)
    
    # ** set parameters
    tol=10.0^(-tol_digits); err_A = err_B = 10000; x=1; # initiate loop variables
    Ãⱼ0 = B̃ᵢ0 = H̃ₘⱼ = H̃ᵣᵢ = w̃ⱼ = ωⱼ = zeros(size(Hₘⱼ,1),1); πᵢⱼ = zeros(size(Hₘⱼ,1),size(Hₘⱼ,1)); # initial guess for model fundamentals
    pos_employment = vec(Hₘⱼ.>0); pos_residence = vec(Hᵣᵢ.>0) # identifying places with firms and residents
    dᵢⱼ = exp.(κ.*τᵢⱼ[findall(pos_employment),findall(pos_residence)]'); H = sum(Hₘⱼ) # by assumption
    while ((err_A<=tol) | (err_B<=tol)) & (x<=200000)
        # equation (12) to recover adjusted wages
        w̃ⱼ[pos_employment] = (Ãⱼ0[pos_employment].^(1/α)).*α.*((1-α)./Qⱼ[pos_employment]).^((1-α)/α); w̃ⱼ[pos_employment] = w̃ⱼ[pos_employment]./geomean(w̃ⱼ[pos_employment]);
        ωⱼ[pos_employment] = w̃ⱼ[pos_employment].^ε; ωⱼ[pos_employment] = ωⱼ[pos_employment]./geomean(ωⱼ[pos_employment]);
        # using equation (4) to estimate commuting probabilities <=> ωⱼ = w̃ⱼᵋ = Eⱼwⱼᵋ ; and B̃ᵢ = BᵢᵋTᵢ
        ϕᵢⱼ = B̃ᵢ0[pos_residence].*(dᵢⱼ.*Qⱼ[pos_residence].^(1-β)).^(-ε).*ωⱼ[pos_employment]';	
        πᵢⱼ[findall(pos_residence),findall(pos_employment)] = ϕᵢⱼ ./ sum(ϕᵢⱼ)
        # using equation (5) to compute predicted workplace and residence employment (remember that L*πᵣᵢ = Hᵣᵢ)
        H̃ᵣᵢ[pos_residence] = sum(πᵢⱼ[findall(pos_residence),findall(pos_employment)], dims=2).*H; 
        H̃ₘⱼ[pos_employment] = sum(πᵢⱼ[findall(pos_residence),findall(pos_employment)], dims=1)'.* H;
        # improving guesses
        Ãⱼ1[pos_employment] = Ãⱼ0[pos_employment] .* (Hₘⱼ[pos_employment]./H̃ₘⱼ[pos_employment]); # rise productivity to rise employment
        B̃ᵢ1[pos_residence] = B̃ᵢ0[pos_residence] .* (Hᵣᵢ[pos_residence]./H̃ᵣᵢ[pos_residence]); # rise amenities to rise residents
        # damping
        Ãⱼ0 = 0.75 .* Ãⱼ0 + 0.25 .* Ãⱼ1; B̃ᵢ0 = 0.75 .* B̃ᵢ0 + 0.25 .* B̃ᵢ1;
        # normalizing for productivity (arbitrary)
        Ãⱼ0[pos_employment] = Ãⱼ0[pos_employment]./geomean(Ãⱼ0[pos_employment]); 
        # normalizing amenities to match data population. 
        "
        The authors measure utility in a unit measure s.t. (Ū/γ)ᵋ/H = 1, where γ = Γ(ε−1/ε) and Γ(·) is the Gamma function (See supplement p. 17).
        Thus, it is implied that ϕ = H, as demonstrated in p. 18 of the supplement. Thus, if the population in the data (H) is greater than the 
        population in the model (ϕ), we increase the amenities to make the city more attractive and attract more residents.
        "
        B̃ᵢ0[pos_residence] = B̃ᵢ0[pos_residence] .* (H./sum(ϕᵢⱼ)).^(1 ./ ε); # if we predict less population than we have in the data, we increase amenities. The ε comes from the B̃ᵢ = BᵢᵋTᵢ definition
        # computing errors
        err_A = round(maximum(abs.(H̃ₘⱼ[pos_employment]-Hₘⱼ[pos_employment])),digits = tol_digits); 
        err_B = round(maximum(abs.(H̃ᵣᵢ[pos_residence]-Hᵣᵢ[pos_residence])),digits = tol_digits);
        x+=1;
        # printing convergence
        println([x, round(err_A/tol, digits=0), round(err_B/tol, digits=0)])
    end
    CMA = sum(ωⱼ' ./ exp.(ν.*τᵢⱼ),dims=2) # equation (29) or section S.3.1.2 
    Ewⱼ = sum((πᵢⱼ ./ sum(πᵢⱼ,dims=2)) .* w̃ⱼ' , 2) # equation (S20)
    return Ãⱼ0, B̃ᵢ0, w̃ⱼ, πᵢⱼ, Ew̃ⱼ, H̃ₘⱼ, H̃ᵣᵢ, CMA
end