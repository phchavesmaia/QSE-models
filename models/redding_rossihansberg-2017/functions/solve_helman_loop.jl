function solveHLwCtyOpen_loop(Lᵢ, wᵢ, aᵢ, dd, noplaces, Iwest, Ieast, H, x)

    # Trade share (we are looking a equation 9.)
    pwmat =  (Lᵢ .* (wᵢ ./ aᵢ) .^ (1-σ)) * ones(1, noplaces)
    nummat = pwmat .* dd
    denom = sum(nummat, dims = 1)
    denommat = ones(noplaces) * denom
    πₙᵢ = nummat ./ denommat

    # Income equals expenditure equilibrium condition (roughly equation 11.)
    income = wᵢ .* Lᵢ
    expend = πₙᵢ * income

    # Domestic trade share
    dom_πₙᵢ = diag(πₙᵢ)

    # Real wage equalization of population and domestic trade share (equation 15.)
    num = ((aᵢ .^ α) .* (H .^ (1-α)) .* (dom_πₙᵢ .^ (-α/(σ-1)))) .^ ((σ-1)/((σ*(1-α))-1))

    # Imposing population immobility across countries (by assumption)
    λₙ = zeros(noplaces)
    Lₑ = zeros(noplaces)
    λₙ[Iwest] = num[Iwest] ./ sum(num[Iwest])
    λₙ[Ieast] = num[Ieast] ./ sum(num[Ieast])
    Lₑ[Iwest] = λₙ[Iwest] .* LLwest
    Lₑ[Ieast] = λₙ[Ieast] .* LLeast

    # Convergence criterion
    income_r = round.(income .* (10^6))
    expend_r = round.(expend .* (10^6))
    Lᵢ_r = round.(Lᵢ .* (10^6))
    Lₑ_r = round.(Lₑ .* (10^6))

    println([x, maximum(abs.(income_r - expend_r)), maximum(abs.(Lₑ_r - Lᵢ_r))])

    # Update loop
    if (income_r == expend_r) & (Lᵢ_r == Lₑ_r)
        display(">>>> Convergence Achieved <<<<")
        x = 10000000
        converge = 1
    else
        wₑ = wᵢ .* (expend ./ income) .^ (1/(σ-1))
        wᵢ = (0.25 .* wₑ) + (0.75 .* wᵢ) # dampening, see https://raw.githack.com/AEM7130/class-repo/master/lecture-notes/04-optimization/04-optimization.html#56
        Lᵢ = (0.25 .* Lₑ) + (0.75 .* Lᵢ) # dampening see https://raw.githack.com/AEM7130/class-repo/master/lecture-notes/04-optimization/04-optimization.html#56

        # Normalization! Choosing geometric mean wage in West as numeraire
        wᵢ[Iwest] = wᵢ[Iwest] ./ geomean(wᵢ[Iwest])
        converge = 0
        x = x + 1
    end

    return x, wᵢ, Lᵢ, πₙᵢ, converge

end