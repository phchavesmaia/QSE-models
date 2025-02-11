function solveHLwCtyOpen_E(fund, d, τⁱ, τᵒ, noplaces)
    " 
    global α, σ, ϕ, LL, LLwest, LLeast
    "
    # Extract location characteristics from fundamentals matrix;
    aᵢ = fund.Aᵢ
    H = fund.Area
    Iwest = fund.Country.==1
    Ieast = fund.Country.==0
    # convergence indicators
    converge = πₙᵢ = 0
    # Initialization based on a symmetric allocation;
    Lᵢ = ones(noplaces) .* (LL/noplaces)
    wᵢ = ones(noplaces)
    # bilateral trade costs
    dd = (d .* (1 .+ τⁱ) .* (1 .+ τᵒ)).^(1-σ)

    # ******************************************************
    # **** START LOOP TO SOLVE FOR WAGES AND POPULATION ****
    # ******************************************************

    x = 1

    while x < 200000
        x, wᵢ, Lᵢ, πₙᵢ, converge = solveHLwCtyOpen_loop(Lᵢ, wᵢ, aᵢ, dd, noplaces, Iwest, Ieast, H, x) 
    end

    # ****************************************************;
    # **** END LOOP TO SOLVE FOR POPULATION AND WAGES ****;
    # ****************************************************;

    return wᵢ, Lᵢ, πₙᵢ, converge

end