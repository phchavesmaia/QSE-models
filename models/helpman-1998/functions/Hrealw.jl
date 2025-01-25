function Hrealw(fund,Lᵢ,πₙᵢ)
    "
    Assumes α, σ, and F as global variables
    Amounts to equation 14.
    "
    πₙₙ = diag(πₙᵢ)

    V̅ = (fund.Aᵢ .^ α) .* (fund.Area .^ (1-α)) .* (πₙₙ .^ (-α/(σ-1))) .* (Lᵢ .^ (-(σ*(1-α)-1)/(σ-1))) 
    V̅ = V̅ ./ (α * ((σ/(σ-1))^α) * ((1/(σ*F))^(α/(1-σ))) * (((1-α)/α)^(1-α)))
    return V̅ 

end
