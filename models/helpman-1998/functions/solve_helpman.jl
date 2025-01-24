function solveHLwCtyOpen_E(fund, d, τⁱ, τᵒ, noplaces)
    " 
    global α, σ, ϕ, LL, LLwest, LLeast
    "
    # set timer
    xtic = time()
    # Extract location characteristics from fundamentals matrix;
    aᵢ = fund.Aᵢ
    H = fund.Area
    #Iwest = fund.Country.==1
    #Ieast = fund.Country.==0
    Iwest = falses(N,N) 
    Ieast = falses(N,N)
    Iwest[:,1:Int(N/2)] .= 1
    Ieast[:,Int(N/2+1):N] .= 1
    Iwest = vec(Iwest)
    Ieast = vec(Ieast)
    # convergence indicators
    converge = πₙᵢ = dom_πₙᵢ = 0
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
        x, wᵢ, Lᵢ, πₙᵢ, dom_πₙᵢ, converge = solveHLwCtyOpen_loop(Lᵢ, wᵢ, aᵢ, dd, noplaces, Iwest, Ieast, H, x) 
    end

    # ****************************************************;
    # **** END LOOP TO SOLVE FOR POPULATION AND WAGES ****;
    # ****************************************************;

    xtic = time() - xtic;

    return wᵢ, Lᵢ, πₙᵢ, dom_πₙᵢ, converge, xtic

end