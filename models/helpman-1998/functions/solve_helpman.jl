function solveHLwCtyOpen_E(fund, d, τⁱ, τᵒ, noplaces)
    " 
    global α, σ, ϕ, LL, LLwest, LLeast
    "
    # set timer
    xtic = time()
    # Extract location characteristics from fundamentals matrix;
    aᵢ = fund.Aᵢ
    H = fund.Area
    Iwest = fund.Country.==1
    Ieast = fund.Country.==0
    # convergence indicators
    convergence = πₙᵢ = dom_πₙᵢ = 0
    # Initialization based on a symmetric allocation;
    Lᵢ = ones(noplaces) .* (LL/noplaces)
    wᵢ = ones(noplaces)
    # bilateral trade costs
    dd = (d .* (1 .+ τⁱ) .* (1 .+ τᵒ)).^(1-σ)

    # ******************************************************
    # **** START LOOP TO SOLVE FOR WAGES AND POPULATION ****
    # ******************************************************

    x = 1
    wind = findall(Iwest .== 1) # returns the cell indices where Iwest is true 
    eind = findall(Ieast .== 1)

    while x < 200000
        # Trade share (we are looking a equation 9.)
        pwmat =  (Lᵢ .* (wᵢ./aᵢ) .^ (1-σ)) * ones(1,noplaces) 
        nummat = pwmat .* dd
        denom = sum(nummat, dims = 1)
        denommat = ones(noplaces) * denom
        πₙᵢ = nummat ./ denommat
        # Income equals expenditure equilibrium condition (roughly equation 11.)
        income = wᵢ .* Lᵢ
        expend = πₙᵢ * income
        # domestic trade share
        dom_πₙᵢ = diag(πₙᵢ)
        # real wage equalization of population and domestric trade share (equation 15.)
        num = ((aᵢ .^ α) .* (H .^ (1-α)) .* (dom_πₙᵢ .^ (-α/(σ-1)))) .^ ((σ-1)/((σ*(1-α))-1))
        # imposing population immobility across countries (by assumption)
        λₙ = Lₑ = zeros(noplaces)
        λₙ[wind] = num[wind] ./ sum(num[wind])
        λₙ[eind] = num[eind] ./ sum(num[eind])
        Lₑ[wind] = λₙ[wind] .* LLwest
        Lₑ[wind] = λₙ[eind] .* LLeast
        
        # Convergence criterion
        income_r = round.(income .* (10^6))
        expend_r = round.(expend .* (10^6))
        Lᵢ_r = round.(Lᵢ .* (10^6))
        Lₑ_r = round.(Lₑ .* (10^6))

        println([x, maximum(abs.(income_r-expend_r)), maximum(abs.(Lₑ_r-Lᵢ_r))])

        # Update loop
        if (income_r == expend_r) & (Lᵢ_r == Lₑ_r)
            display(">>>> Convergence Achieved <<<<")
            x = 10000000
            converge = 1;
        else
            wₑ = wᵢ .* (expend./income) .^ (1/(σ-1)) # Idk where it comes from
            wᵢ = α .* wₑ + (1-α) .* wᵢ # Idk, probably wouldn't matter
            Lᵢ = α .* Lₑ + (1-α) .* Lᵢ # Idk, probably wouldn't matter

            # Normalization! Choosing geometric mean wage in West as numeraire
            wᵢ = wᵢ ./ geomean(wᵢ[wind])
            converge = 0
            x=x+1;
        end
    end

    # ****************************************************;
    # **** END LOOP TO SOLVE FOR POPULATION AND WAGES ****;
    # ****************************************************;

    xtic = time() - xtic;

    return wᵢ, Lᵢ, πₙᵢ, dom_πₙᵢ, converge, xtic

end