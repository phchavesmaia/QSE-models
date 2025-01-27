function Hpindex(fund,Lᵢ,wᵢ,πₙᵢ) 
    "
    Assumes α, σ, LL, and F as global variables
    This function amounts to equation 8.
    "
    aᵢ=fund.Aᵢ ; πₙₙ = diag(πₙᵢ);  

    return (σ/(σ-1)).*((Lᵢ./((σ*F).*πₙₙ)).^(1/(1-σ))).*(wᵢ./aᵢ)

end