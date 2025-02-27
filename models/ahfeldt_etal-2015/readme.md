# Ahfeldt et al. (2015)

The goal of this folder is to replicate the canonical QSE model from Ahlfeldt et al. (2015) using Julia, based on the Matlab Toolkit by Ahlfeldt (2024). To download the .mat files, please refer to Ahlfeldt (2024).

Since `GeoIO.jl' doesn't support EPSG:3068, I reproject the data to EPSG:3035 using QGIS to ensure a seamless workflow. These altered data are the shapefiles with "1" at the end of their names.

It follows a list of what remains to be implemented:
- [ ] Test whether the results of the calibration function return identical results to Ahlfeldt (2024).
    - $\tilde{A}_j$ calibration is clearly wrong, suggesting estimation problems wrt $\omega_j$. 
- [ ] Implement the counterfactual exercises with **exogenous** fundamentals.
- [ ] Implement the counterfactual exercises with **endogenous** fundamentals

## References

Ahlfeldt, G. M. (2024). Toolkit for The Economics of Density: Evidence from the Berlin Wall. https://github.com/Ahlfeldt/ARSW2015-toolkit

Ahlfeldt, G. M., Redding, S. J., Sturm, D. M., & Wolf, N. (2015). The Economics of Density: Evidence from the Berlin Wall. Econometrica, 83(6), 2127-2189. https://doi.org/10.3982/ECTA10876