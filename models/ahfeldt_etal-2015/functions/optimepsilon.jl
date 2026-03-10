function get_ω(Hₘⱼ,Hᵣᵢ,τᵢⱼ,Qⱼ; tol_digits=6, x_max = 500, ε = 1)
    "
    This function solves for TRANSFORMED wages (ωⱼ) for given values of
    workplace employment, residence employment, and bilateral commuting
    costs. If ε is not provided, it will be set to 1 so that you can 
    compute the initial guess of ω without having to compute ε first.
    "
    # initiate main loop and output variables
    n_places = size(Hₘⱼ,1);
    ωⱼ  = zeros(n_places,1); Ĥₘⱼ = zeros(n_places,1);
    pos_employment = vec(Hₘⱼ.>0); pos_residence = vec(Hᵣᵢ.>0) ; # identifying places with firms and residents
    x=1; err = 10000; tol = 10.0^(-tol_digits); # defining loop variables

    # now, I ONLY care for places that are being used 
    τᵢⱼ = τᵢⱼ[findall(pos_residence),findall(pos_employment)] ; 
    πᵢⱼi = zeros(size(τᵢⱼ)) ; 
    Hᵣᵢ = Hᵣᵢ[pos_residence]; 
    Hₘⱼ = Hₘⱼ[pos_employment]; n_workplaces = size(Hₘⱼ,1);

    # initiating some more variables for numerical efficiency
    Ĥₘⱼ0 = zeros(n_workplaces,1); 
    dᵢⱼε = @. exp(ν * τᵢⱼ); # pre-computing the, by assumption, iceberg commuting cost times ε for numerical efficiency
    
    # initial guess on transformed wages ωⱼ 
    ωⱼ0 = ωⱼ[pos_employment]; 
    w̃ⱼ0 = @. (((1-α)/Qⱼ[pos_employment])^((1-α)/α))*α # initial guess on ADJUSTED wages following Equation (12) which combines first-order condition and zero-profit conditions, after setting Aⱼ = 1.
    @. ωⱼ0 = w̃ⱼ0 ^ ε; # initial guess on TRANSFORMED wages using that ω = w̃^ε.
    ωⱼ1 = zeros(n_workplaces,1);

    # announcing the function
    println(">>>> Calibrating ω <<<<")
    while (err >= tol) && (x <= x_max)
        # Compute conditional commuting probabilities (equation 4)
        @. πᵢⱼi = (ωⱼ0' / dᵢⱼε) / $sum(ωⱼ0' ./ dᵢⱼε, dims=2) ;
        # Compute predicted workplace employment (equation 7 or, more explicitly, equation 26 and S.44)
        @. Ĥₘⱼ0 = $sum(πᵢⱼi * Hᵣᵢ, dims=1)' ;
        # Compute Employment Gap and Check Convergence
        err = @. $round($maximum(abs(Ĥₘⱼ0 - Hₘⱼ)),digits = tol_digits) ;
        # Update ω guess
        @. ωⱼ1 = ωⱼ0 * (Hₘⱼ / Ĥₘⱼ0) ;
        # Apply damping to improve stability (I will follow ARSW and use a 0.5 damping factor, even if 0.75/0.25 should be safer)
        @. ωⱼ0 = 0.5 * ωⱼ0 + 0.5 * ωⱼ1 ;
        # Normalize wages to ensure geomean(ωⱼ) = 1
        @. ωⱼ0 = ωⱼ0 ./ $geomean(ωⱼ0);
        # Print convergence rate
        println([x, trunc(err / tol, digits=0)])
        x += 1;
    end
    if x==x_max
        error("Convergence not achieved for adjusted wages (ω)")
    end
    
    @. ωⱼ[pos_employment] = ωⱼ0
    @. Ĥₘⱼ[pos_employment] = Ĥₘⱼ0
    println(">>>> Wage System Converged <<<<")

    return ωⱼ, Ĥₘⱼ
end

function payroll_aggregator(su)
    "
    The su (spatial unit) variable is a map between smaller
    spatial units, such as blocks, to larger spatial
    units (lsu), such as districts. Each line of `spatial_unit`
    regards a specific su in accordance to the row number, 
    whereas its values correspond to the lsu. 
    "
    lsu = unique(su);
    n_su = length(su);

    # indexing the lsu
    lsu_map = Dict(id => i for (i,id) in enumerate(lsu));

    # building sparse matrix A where A[lsu_index,su_index] = 1, i.e., 
    # it indicates 1 if a su is part of a lsu. Naturally, A is n_lsu x n_su.
    I = vec([lsu_map[id] for id in su]); # translates su values (lsu code) to index (lsu_map values)
    J = 1:n_su;
    V = ones(n_su);
    S = sparse(I, J, V); 
    return S
end

function get_fϵ(u,p)
    "
    Defining the objective function of the minimization problem
    to find ε.
    "
    
    # *************************
    # *** Unpack parameters ***
    # *************************
    ε = u[1];
    S, Hₘⱼ, ωⱼ, Vlwⱼ = p;

    # *******************
    # ****** Wages ******
    # *******************

    # compute payroll at the block level
    w̃ⱼ = @. ωⱼ ^ (1/ε);
    @. w̃ⱼ[w̃ⱼ.>0] = w̃ⱼ[w̃ⱼ.>0] / $geomean(w̃ⱼ[w̃ⱼ.>0]); # normalizing after the change
    payroll_su = @. w̃ⱼ * Hₘⱼ;

    # aggregating payroll and labor to lsu levels
    payroll_lsu = S * payroll_su;
    labor_lsu = S * Hₘⱼ;

    # getting wages at the lsu level
    w̃ⱼ_lsu = @. payroll_lsu / labor_lsu;
    lw̃ⱼ_lsu = log.(w̃ⱼ_lsu);                                                  
    @. lw̃ⱼ_lsu = lw̃ⱼ_lsu - $mean(lw̃ⱼ_lsu); # demean                                                
    Vlw̃ⱼ_lsu = var(lw̃ⱼ_lsu);

    # *******************************
    # ****** Moment Conditions ******
    # *******************************

    ftD = Vlw̃ⱼ_lsu - Vlwⱼ; # error
    ftt = ftD^2 * 10.0^6; # square error (multiplied for numerical consistency), equivalent to equation S.64
    "
    Observe that ftt (equivalent to equation 35 or S.64), which should be 0, can be read as:
    E[(1/ε)²⋅log(ω)² - σₗₙ₍w₎²] = 0
    Thus, we can use this moment condition to identify ε as it is the only unkown in the equation.
    Notice further that E[(1/ε)²⋅log(ω)²] is the variance of transformed wages 
    since ω has a mean of 1 and, hence, ln(ω)=0.
    "
    return ftt
end

function get_ε(Vlwⱼ,Hₘⱼ,Hᵣᵢ,τᵢⱼ,Qⱼ; su=vec(1:length(Qⱼ)), tol_digits = 6, ε0=4, maxiter=1000)
    
    # *****************************************
    # ******* Computing ajusted wages ω *******
    # *****************************************

    ωⱼ, Ĥₘⱼ = get_ω(Hₘⱼ,Hᵣᵢ,τᵢⱼ,Qⱼ,tol_digits=tol_digits);    

    # ************************************************
    # ******* Defining optimazation parameters *******
    # ************************************************

    S = payroll_aggregator(su);
    p = (S, Hₘⱼ, ωⱼ, Vlwⱼ);
    u0 = [ε0];

    # ***************************
    # ******* Computing ε *******
    # ***************************
    
    # Define the Problem 
    # --- OptimizationProblem(function, initial_guess, parameters; lb, ub)
    prob = OptimizationProblem(get_fϵ, u0, p; lb = [2.0], ub = [24.0]);
    # Solve the problem 
    # --- We use the BOBYQA algorithm, differently from the original implementation that used the Generalized Pattern Search (GPS) algorithm.
    println(">>>> Calibrating ε <<<<")
    sol = solve(prob, NLopt.LN_BOBYQA(), reltol = 10.0^(-tol_digits));
    println(
    """
    objective value       : $(sol.objective)
    solution (ε)          : $(sol.u[1])
    solution status       : $(sol.retcode)
    # function evaluation : $(sol.stats.fevals)
    """
    )
    return sol.u[1], Ĥₘⱼ, ωⱼ
end