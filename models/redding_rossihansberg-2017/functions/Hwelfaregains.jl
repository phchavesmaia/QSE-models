function Hwelfaregains(πₙᵢ,πₙᵢᶜ,Lᵢ,Lᵢᶜ)
    "
    α and σ are assumed to be global variables // NOT WORKING!!
    "
    # defining domestic trade shares
    πₙₙ = diag(πₙᵢ); πₙₙᶜ = diag(πₙᵢᶜ);

    V̅̂ = ((πₙₙ./πₙₙᶜ).^(α/σ-1)).*((Lᵢ./Lᵢᶜ).^((σ*(1-α)-1)/(σ-1))) # this is equation 21.
    
    return V̅̂
end