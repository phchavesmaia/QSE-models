function get_ω_and_ε(Vlwⱼ,Hₘⱼ,Hᵣᵢ,τᵢⱼ; tol_digits = 6, ε0=4, maxiter=1000)
    
    # *****************************************
    # ******* Computing ajusted wages ω *******
    # *****************************************

    function get_ω(Hₘⱼ,Hᵣᵢ,τᵢⱼ ; tol_digits=tol_digits, x_max = 500)
        ωⱼ = zeros(size(Hₘⱼ,1),1); x=1; 
        pos_employment = vec(Hₘⱼ.>0); pos_residence = vec(Hᵣᵢ.>0) # identifying places with firms and residents
        err = 10000; tol = 10.0^(-tol_digits); # defining loop variables

        # ** now, I ONLY care for places that are being used ** 
        τᵢⱼ = τᵢⱼ[findall(pos_employment),findall(pos_residence)]' ; Hᵣᵢ = Hᵣᵢ[pos_residence]; Hₘⱼ = Hₘⱼ[pos_employment];
        # ** initial guess on ω 
        ωⱼ0 = ωⱼ[pos_employment]; # Again, on wages, I only care for places with firms i.e. that pay wages (workplace employment)
        ωⱼ0=(((1-α)./Qⱼ[pos_employment]).^((1-α)/α)).*α; # Equation (12) which combines first-order condition and zero-profit conditions.

        # announcing the function
        println(">>>> Calibrating ω <<<<")
        while (err>=tol) & (x<=x_max) 
            # predicting workplace employment (eq. 26 or appendix S.44)
            Ĥₘⱼ = sum((ωⱼ0' ./ exp.(ν.*τᵢⱼ).* Hᵣᵢ)./sum(ωⱼ0' ./ exp.(ν.*τᵢⱼ), dims=2), dims=1)';
            # calculate error
            err = round(maximum(abs.(Ĥₘⱼ-Hₘⱼ)),digits = tol_digits)
            # update ω guess
            ωⱼ1 = ωⱼ0 .* (Hₘⱼ./Ĥₘⱼ)
            # damping 
            ωⱼ0 = 0.75 .* ωⱼ0 + 0.25 .* ωⱼ1 
            # normalizing wages to geomean(ωⱼ) = 1
            ωⱼ0 = ωⱼ0./geomean(ωⱼ0)
            # print convergence
            #println([x, round(err/tol, digits=0)])
            x+=1;
        end
        if x==x_max
            error("Convergence not achieved for adjusted wages (ω)")
        end
        ωⱼ[pos_employment] = ωⱼ0
        println(">>>> Wage System Converged <<<<")
        return ωⱼ
    end

    # **********************************************************
    # ******* Computing value of objective function f(ε) *******
    # **********************************************************

    ωⱼ = get_ω(Hₘⱼ,Hᵣᵢ,τᵢⱼ);

    function get_fϵ(ε)
        # *******************
        # ****** Wages ******
        # *******************
        w̃ⱼ = ωⱼ .^ (1/ε);
        w̃ⱼ[w̃ⱼ.>0] = w̃ⱼ[w̃ⱼ.>0]./geomean(w̃ⱼ[w̃ⱼ.>0]) # normalizing after the change
        payroll = w̃ⱼ.* Hₘⱼ # generate the payrolls
        # aggregating payroll to Bezirke level
        df = DataFrame([payroll block_bzk], :auto)
        grouped_df = combine(groupby(df, :x2), :x1 => mean => :mean_value, :x1 => length => :count)
        payroll_bzk = grouped_df.mean_value .* grouped_df.count
        # aggregating jobs to Bezirke level
        df = DataFrame([Hₘⱼ block_bzk], :auto)
        grouped_df = combine(groupby(df, :x2), :x1 => mean => :mean_value, :x1 => length => :count)
        labor_bzk = grouped_df.mean_value .* grouped_df.count
        # getting wages at the Bezirke level
        w̃ⱼbzk = payroll_bzk./labor_bzk
        lw̃ⱼbzk=log.(w̃ⱼbzk)                                                  
        lw̃ⱼbzk=lw̃ⱼbzk.-mean(lw̃ⱼbzk)                                                
        Vlw̃ⱼbzk=var(lw̃ⱼbzk)
        # *******************************
        # ****** Moment Conditions ******
        # *******************************
        ftD = Vlw̃ⱼbzk - Vlwⱼ; # error
        ftt = ftD^2 .* 10.0^6; # square error (multiplied for numerical consistency), equivalent to equation S.64
        "
        Observe that ftt (equivalent to equation 35 or S.64), which should be 0, can be read as:
        E[(1/ε)²⋅log(ω)² - σₗₙ₍w₎²] = 0
        Thus, we can use this moment condition to identify ε as it is the only unkown in the equation.
        Notice further that E[(1/ε)²⋅log(ω)²] is the variance of transformed wages 
        since ω has a mean of 1 and, hence, ln(ω)=0.
        "
        return ftt
    end

    # ***************************
    # ******* Computing ε *******
    # ***************************

    # defining the optimization algorithm. 
    "
    We use the BOBYQA algorithm, differently from the original implementation that used the Generalized Pattern Search (GPS) algorithm.
    "
    println(">>>> Calibrating ε <<<<")
    opt = Opt(:LN_BOBYQA, 1)
    lower_bounds!(opt, [2])
    upper_bounds!(opt, [24])
    xtol_rel!(opt, 10.0^(-tol_digits))
    min_objective!(opt, (x,grad) -> get_fϵ(x[1]))  # so that we can ignore the package gradient requirement
    maxeval!(opt, maxiter)
    (minf, minx, ret) = optimize(opt, [ε0])
    num_evals = NLopt.numevals(opt)
    println(
    """
    objective value       : $minf
    solution (ε)          : $minx
    solution status       : $ret
    # function evaluation : $num_evals
    """
    )
    return ωⱼ, minx[1]
    
end