function Hlandprice(fund,Lᵢ,wᵢ)
    "
    Assumes α, σ, LL, LLwest, and LLeast as global variables
    Amounts to equation 12.
    "
    Hᵢ = fund.Area; 
    
    return ((1-α)/α).* (wᵢ.*Lᵢ ./ Hᵢ)

end