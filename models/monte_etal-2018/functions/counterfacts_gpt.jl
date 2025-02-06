function counterFacts_gpt(wₙ,v̄ₙ,Lₙ,Rₙ,πₙᵢ,λₙᵢ; Âᵢ=ones(size(wₙ)), B̂ₙᵢ=ones(size(πₙᵢ)), κ̂ₙᵢ=ones(size(πₙᵢ)), d̂ₙᵢ=ones(size(πₙᵢ)), tol_digits = 4)
    "
    This function executes a series of counterfactuals in which we 
    simulate the effects of a hypothetical border that follows the 
    border between formerly separated West Germany and East Germany 

    The necessary arguments are the model variables, whereas the kwds
    denote the exact hat algebra counterfactual exercise. If a kwd is
    not determined, it is set to ones (that primitive does not change).

    wₙ, v̄ₙ, Lₙ, and Rₙ are Nx1 vectors (or 1 column matrices)
    πₙᵢ and λₙᵢ are NxN matrices
    "

    # Define main loop parameters
    err_ŵₙ = err_λ̂ₙᵢ = 1e4
    tol = 10.0^(-tol_digits)
    x = 1

    # Initialize guesses of changes
    ŵₙ0 = ones(size(wₙ))
    λ̂ₙᵢ0 = ones(size(λₙᵢ))

    # Define auxiliary functions
    num_v̄̂ₙ(λₙᵢ, ŵₙ, B̂ₙᵢ, κ̂ₙᵢ) = λₙᵢ .* B̂ₙᵢ .* (ŵₙ' ./ κ̂ₙᵢ) .^ ε
    num_π̂ₙᵢ(L̂ₙ, ŵₙ, Âᵢ, d̂ₙᵢ) = L̂ₙ' .^ (1 - (1 - σ) * ν) .* (d̂ₙᵢ .* ŵₙ' ./ Âᵢ') .^ (1 - σ)
    num_λ̂ₙᵢ1(B̂ₙᵢ, P̂ₙ, P̂ₕₙ, ŵₙ, κ̂ₙᵢ) = B̂ₙᵢ .* (P̂ₙ .^ α .* P̂ₕₙ .^ (1 - α)) .^ (-ε) .* (ŵₙ' ./ κ̂ₙᵢ) .^ ε

    # Define inside-loop variables so that they "survive" loop completion
    local v̄̂ₙ, L̂ₙ, R̂ₙ, P̂ₕₙ, π̂ₙᵢ, P̂ₙ

    while (err_ŵₙ >= tol) || (err_λ̂ₙᵢ >= tol)
        "
        All of these updating equations can be found in SW appendix A.2. or MRRH appendix B
        "
        # Calculate v̄̂ₙ
        v̄̂ₙ = (1 ./ v̄ₙ) .* sum(num_v̄̂ₙ(λₙᵢ, ŵₙ0, B̂ₙᵢ, κ̂ₙᵢ) ./ sum(num_v̄̂ₙ(λₙᵢ, ŵₙ0, B̂ₙᵢ, κ̂ₙᵢ), dims=1) .* (ŵₙ0 .* wₙ), dims=1)

        # Calculate L̂ₙ (employment)
        L̂ₙ = sum(Lₙ) ./ Lₙ .* (sum(λₙᵢ .* λ̂ₙᵢ0, dims=2))  # Corrected sum alignment

        # Calculate R̂ₙ (residence)
        R̂ₙ = sum(Lₙ) ./ Rₙ .* sum(λₙᵢ .* λ̂ₙᵢ0, dims=1)  # Ensuring consistency with MATLAB

        # Calculate P̂ₕₙ (house prices)
        P̂ₕₙ = (v̄̂ₙ .* R̂ₙ) .^ (1 / (1 + δ))

        # Calculate π̂ₙᵢ (trade shares)
        π̂ₙᵢ = num_π̂ₙᵢ(L̂ₙ, ŵₙ0, Âᵢ, d̂ₙᵢ)' ./ sum(πₙᵢ .* num_π̂ₙᵢ(L̂ₙ, ŵₙ0, Âᵢ, d̂ₙᵢ), dims=1)

        # Calculate P̂ₙ (price index)
        P̂ₙ = (L̂ₙ .^ (1 - (1 - σ) * ν) ./ diag(π̂ₙᵢ)) .^ (1 / (1 - σ)) .* (diag(d̂ₙᵢ) .* ŵₙ0 ./ Âᵢ)

        # Calculate loop variables ŵₙ1 and λ̂ₙᵢ1
        ŵₙ1 = (1 ./ (wₙ .* Lₙ .* L̂ₙ)) .* sum(πₙᵢ .* π̂ₙᵢ .* v̄ₙ .* v̄̂ₙ .* Rₙ .* R̂ₙ, dims=2)
        ŵₙ1 = ŵₙ1 ./ mean(ŵₙ1 .* wₙ)  # Normalization step to match MATLAB
        λ̂ₙᵢ1 = num_λ̂ₙᵢ1(B̂ₙᵢ, P̂ₙ, P̂ₕₙ, ŵₙ0, κ̂ₙᵢ) ./ sum(λₙᵢ .* num_λ̂ₙᵢ1(B̂ₙᵢ, P̂ₙ, P̂ₕₙ, ŵₙ0, κ̂ₙᵢ), dims=1)

        # Calculate errors
        err_ŵₙ = round(maximum(abs.(ŵₙ1 - ŵₙ0)), digits=tol_digits)
        err_λ̂ₙᵢ = round(maximum(abs.(λ̂ₙᵢ1 - λ̂ₙᵢ0)), digits=tol_digits)

        # Update loop variables
        ŵₙ0 = 0.25 * ŵₙ1 + 0.75 * ŵₙ0
        λ̂ₙᵢ0 = 0.25 * λ̂ₙᵢ1 + 0.75 * λ̂ₙᵢ0

        # Print convergence rate
        println([x, err_ŵₙ, err_λ̂ₙᵢ])
        x += 1
    end

    # Calculate welfare changes
    Ū̂ = ŵₙ0' .* (B̂ₙᵢ ./ λ̂ₙᵢ0) .^ (1 / ε) ./ (κ̂ₙᵢ .* P̂ₙ .^ α .* P̂ₕₙ .^ (1 - α))
    percentageChange = round((Ū̂[1, 1] - 1) * 100, digits=2)
    println("...Change in welfare is " * string(percentageChange) * "%")

    return ŵₙ0, v̄̂ₙ, P̂ₕₙ, π̂ₙᵢ, λ̂ₙᵢ0, P̂ₙ, R̂ₙ, L̂ₙ, Ū̂
end
 