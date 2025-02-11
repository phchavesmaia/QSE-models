function get_bi_lamba(wₙ,τₙᵢ,Lᵢ,Rₙ,L̄; tol_digits=4)
    " 
    τₙᵢ denotes commuting time
    "

    J = size(wₙ,1); Tᵢ0 = ones(J,1) # initial guess
    err_b = err_t = 10000; x = 1; tol = 10.0^(-tol_digits); κₙᵢ = τₙᵢ.^ μ # loop variables
    local λₙᵢn # initiating variables that I want to keep outside the loop

    while (err_t>=tol) & (x<=200000)

        # getting conditional commuting probabilities from eq. 12 of MRRH, where I assume Bₙᵢ = Tₙ * Eᵢ (residential and workplace amenities, following the ARSW notation)
        λₙᵢn = (Tᵢ0' .* (wₙ' ./ κₙᵢ) .^ ε) ./ sum(Tᵢ0' .* (wₙ' ./ κₙᵢ) .^ ε, dims= 2)
        Lᵢ¹ =  sum(λₙᵢn .* Rₙ, dims=1)' # using eq 13 of MRRH to get labor that rationalizes λₙᵢi
        # getting Tᵢ1
        Tᵢ1 = Tᵢ0 .* (Lᵢ./Lᵢ¹) # increase workplace amenity to attract more workers if observed employment exceeds predicted employment
        
        # get loop error
        err_t = round(maximum(abs.(Tᵢ1-Tᵢ0)),digits=tol_digits)

        # damping
        Tᵢ0 = 0.75.*Tᵢ0 + 0.25.*Tᵢ1

        # print convergence rate
        println([x, round(err_t/tol)])
        x+=1
    end

    # Compute residential and workplace choice probabilities (equation 11)
    λₙᴿ = Rₙ/sum(Rₙ)

    # Compute unconditional location choice probabilities (sorta eq. 12 again)
    λₙᵢ = λₙᵢn .* λₙᴿ

    # Compute commuting flows
    Lₙᵢ = λₙᵢ .* L̄

    return Tᵢ0, λₙᵢn, λₙᵢ, Lₙᵢ
end