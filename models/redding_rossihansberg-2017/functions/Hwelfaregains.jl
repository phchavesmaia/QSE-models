function Hwelfaregains(πₙᵢ,πₙᵢᶜ,L̂ᵢ)
    "
    α and σ are assumed to be global variables
    "
    # defining domestic trade shares
    πₙₙ = diag(πₙᵢ); πₙₙᶜ = diag(πₙᵢᶜ); π̂ₙₙ = πₙₙᶜ ./ πₙₙ;

    V̅̂ = ((1 ./ π̂ₙₙ).^(α/(σ-1))).*((1 ./ L̂ᵢ).^((σ*(1-α)-1)/(σ-1))) # this is equation 21.
    
    return V̅̂
end