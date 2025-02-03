function solveProductTrade(Lₙ, Rₙ, wₙ, v̄ₙ, dₙᵢ)

    "
    This function recovers fundamental productivities that satisfy the  
    equals expenditure condition and uses them to comute trade   
    trade shares and the tradable goods price index  
    
    This is fundamentally based on Eqs. (10) and (12) of SW
    "
    Aᵢ0 = ones(J);
    # defining numerator of eq (10):   
    numerator(Lₙ,wₙ,dₙᵢ,Aᵢ) = (Lₙ'.^(1-(1-σ)*ν)) .* ((wₙ' .* dₙᵢ ./ Aᵢ') .^ (1-σ)) # observe that labor, wages and productivity change across columns!
    # defining main loop
    err=10000; tol = 1e-6; x = 1; 
    while (err>=tol) & (x<=200000)

        # define trade shares
        πₙᵢ = numerator(Lₙ,wₙ,dₙᵢ,Aᵢ0) ./ sum(numerator(Lₙ,wₙ,dₙᵢ,Aᵢ0), dims=2)
        
        # defining income and expenditure, i.e., equation (12) components:
        income = wₙ .* Lₙ; expenditure = sum(πₙᵢ .* v̄ₙ .* Rₙ, dims=1)';
        
        # equation (12) -- equilibrium condition
        err = round.(maximum(abs.(income - expenditure)),digits=6);

        # Update Aᵢ guess
        Aᵢ1 = Aᵢ0 .* (income ./ expenditure) # if income>expenditure => workers need to be more productive to afford their purchases => higher Aᵢ 
        
        # Damping to update the value and update iteration number
        Aᵢ0 = (0.25 .* Aᵢ1) + (0.75 .* Aᵢ0); x+=1;
        
        #Productivity is identified up to a constant. Therefore we normalize to by the mean
        Aᵢ0 = Aᵢ0 ./ mean(Aᵢ0) 

        # report convergence ratio
        #println([x, err])

    end

    # converged?
    if x==200000
        error("Convergence not achieved")
    end

    # get prices
    Pₙ = (σ/(1-σ)) .* (((Lₙ.^(1-(1-σ)*ν))./(σ.*f.*diag(πₙᵢ))) ) .^ (1/(1-σ)) .* (wₙ .* diag(dₙᵢ) ./ Aᵢ0)

    return Aᵢ0, πₙᵢ, Pₙ

end