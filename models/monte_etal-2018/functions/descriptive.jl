function descriptive_analysis(Lₙ, Areaₙ, Rₙ, distₙᵢ, comMat, Aᵢ, baseline, λₙᵢ)    
    
    # Compute densities
    lEmpDensity_n = log.(Lₙ) .- log.(Areaₙ)
    lPopDensity_n = log.(Rₙ) .- log.(Areaₙ)
    lCommImport_n = log.(Lₙ) .- log.(Rₙ)


    # Map some input data

    mapit("./data/shapes/VG250_KRS_clean_final.shp", lEmpDensity_n, "Log employment density", path_to="./figures/MAP_EmpDensity.png")
    mapit("./data/shapes/VG250_KRS_clean_final.shp", lPopDensity_n, "Log population density", path_to="./figures/MAP_PopDensity.png")
    mapit("./data/shapes/VG250_KRS_clean_final.shp", lCommImport_n, "Log commuting import", path_to="./figures/MAP_CommImport.png")
    mapit("./data/shapes/VG250_KRS_clean_final.shp", Pₕₙ, "House price", path_to="./figures/MAP_HousePrice.png", label_legend = "Rent index")
    mapit("./data/shapes/VG250_KRS_clean_final.shp", wₙ, "Wage", path_to="./figures/MAP_Wage.png", label_legend = "Wage index")

    # Check distances
    avDist = mean(distₙᵢ, dims=2)
    mapit("./data/shapes/VG250_KRS_clean_final.shp", avDist, "Average distance to other counties", path_to="./figures/MAP_AvDist.png", label_legend = "Distance index")

    # Correlate commuting probabilities with distance
    scatter(log.(distₙᵢ), log.(comMat), xlabel="Distance (log)", ylabel="Commuting probability (log)", title="Commuting probability, log-log model", grid=true, legend=false)
    savefig("./figures/Scatter_Lambda_log_log.png")
    scatter(distₙᵢ, log.(comMat), xlabel="Distance (km)", ylabel="Commuting probability (log)", title="Commuting probability, log-lin model", grid=true, legend=false)
    savefig("./figures/Scatter_Lambda_log_lin.png")

    # Fixed effects regressions
    comMat_vec = vec(comMat)
    dist_mat_vec = vec(distₙᵢ)

    # Ensuring all entries are positive before taking logs
    positive_idx = comMat_vec .> 0
    log_comMat_vec = log.(comMat_vec[positive_idx])
    log_dist_mat_vec = log.(dist_mat_vec[positive_idx])
    dist_mat_vec = dist_mat_vec[positive_idx]
    log_λₙᵢ_vec = log.(vec(λₙᵢ)[positive_idx])
    log_time_vec = log.(vec(baseline)[positive_idx])

    # plotting scatterplot of positive coomutes data
    scatter(dist_mat_vec, log_comMat_vec, xlabel="Distance (Km)", ylabel="Commuting probability (log)", title="Commuting probability, log-log model", label=false, grid=true)
    savefig("./figures/Scatter_Lambda_log_lin_Positive_Commute.png")

    # Generate dummy variables
    n = size(distₙᵢ, 1)
    df = DataFrame(reduce(hcat,collect.(Tuple.(findall(reshape(positive_idx, n, n)))))', :auto)
    rename!(df,[:row_indices, :col_indices])
    df[!,"log_comMat_vec"] = log_comMat_vec
    df[!,"log_dist_mat_vec"] = log_dist_mat_vec
    df[!,"log_λₙᵢ_vec"] = log_λₙᵢ_vec
    df[!,"log_time_vec"] = log_time_vec
    df[!,"dist_mat_vec"] = dist_mat_vec

    model = reg(df, @formula(log_λₙᵢ_vec ~ log_dist_mat_vec  + fe(row_indices) + fe(col_indices))) # row=origin; col=destination
    slope = coef(model)[1]
    println("The slope of the regression is: $slope")

    # For log-lin model
    model = reg(df, @formula(log_λₙᵢ_vec ~ dist_mat_vec  + fe(row_indices) + fe(col_indices))) 
    Comm_slope = coef(model)[1]
    println("The slope of the regression is: $Comm_slope")

    # Compute commuter market access
    ωₙ = wₙ .^ ε
    CommWeight_ni = distₙᵢ .^ (-μ .* ε)
    CMA_n = CommWeight_ni * ωₙ
    mapit("./data/shapes/VG250_KRS_clean_final.shp", log.(CMA_n), "Log commuting market access", path_to="./figures/MAP_CMA.png")
    EmpPot_n = CommWeight_ni * Lₙ
    mapit("./data/shapes/VG250_KRS_clean_final.shp", log.(EmpPot_n), "Log employment potential", path_to="./figures/EmpPot.png")

    # Map productivities
    mapit("./data/shapes/VG250_KRS_clean_final.shp", log.(Aᵢ), "Log productivity", path_to="./figures/MAP_A.png" )
    scatter(log.(Aᵢ), log.(wₙ), xlabel="Log productivity A", ylabel="Log wage w", title="Productivity vs. wage", grid=true, legend=false)
    savefig("./figures/Scatter_A_w.png")

    # Map trade shares and price index
    mapit("./data/shapes/VG250_KRS_clean_final.shp", diag(πₙᵢ), "Own trade share", path_to="./figures/MAP_pi_nn.png", label_legend = "Percentage points")
    mapit("./data/shapes/VG250_KRS_clean_final.shp", Pₙ, "Tradables price index", path_to="./figures/MAP_P_n.png", label_legend = "Price index")

    return println("<<<<<<<<<<<<<<< Descriptive analysis completed >>>>>>>>>>>>>>>")

end