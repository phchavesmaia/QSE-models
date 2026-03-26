# Ahfeldt et al. (2015)

The goal of this folder is to replicate the canonical QSE model from Ahlfeldt et al. (2015) using Julia, based on the Matlab Toolkit by Ahlfeldt (2024). To download the .mat files, please refer to Ahlfeldt (2024). Importantly, I renamed `prepdata_big_TD.mat` to `prepdata_big_TD06.mat` to make it easier for the helping import functions.

Since `GeoIO.jl' doesn't support EPSG:3068, I reproject the data to EPSG:3035 using QGIS to ensure a seamless workflow. These altered data are the shapefiles with "1" at the end of their names.

It follows a list of what remains to be implemented:
- [X] Test whether the results of the calibration function (model inversion) return identical results to Ahlfeldt (2024).
    - [X] Sequential algorithim;
    - [X] Simultaneous algorithm;
    - [X] Both algorithms return identical results.
- [X] Implement the counterfactual exercises with **exogenous** fundamentals.
    - [X] Solve model for exogenous fundamentals;
    - [X] Implement counterfactuals.
- [X] Implement the counterfactual exercises with **endogenous** fundamentals
    - [X] Solve model for endogenous fundamentals;
    - [X] Implement counterfactuals.
        - [X] Closed-city
        - [X] Open-city

## Notation
In this replication, I opted for a different approach compared to the one used for the other models. Here, instead of writing the point-wise operations as:
```julia
Y = X .* Z .+ K
```
I wrote it as:
```julia
Y = @. X * Z + K
```
That is, I used the `@ .` broadcasting macro as much as possible. This is very important for the betterment efficiency and memory allocation. When I am defining a new variable, I use the macro after the `=` sign, as in the example above. When changing the values of a pre-defined variable, I use the macro before the whole line, as in:
```julia
@. Y = X * Z + K
```
where `Y` was previously defined in the code. This is correct way to guarantee maximum benefits. 

Importantly, the macro tries to broadcast and fuse all operations it contemplates, which may lead to problems such as broadcasting the `mean()` function. To avoid it, I flair these functions with symbol `$` to scape the broadcasting, as in the following example:
```julia
@. Y = X * Z + $mean(K)
``` 

## Structure
In this repository, I slightly reorganized the folder structure. Instead of creating a `functions` folder to be directly included in the main file, I introduced a `modules` folder. This change was intentional, as I wanted to learn how to properly build and integrate modules within a Julia workflow.

This approach makes the project more modular overall. Each module explicitly declares its own dependencies, which keeps the main file clean and makes it easier to reuse specific pieces of code in other projects.

The trade-off is that the code becomes more “Julian,” in the sense that it can be a bit harder to navigate at first. If you prefer a simpler structure without modules, I recommend checking out a commit from March 11, 2026, or earlier. 

## References

Ahlfeldt, G. M. (2024). Toolkit for The Economics of Density: Evidence from the Berlin Wall. https://github.com/Ahlfeldt/ARSW2015-toolkit

Ahlfeldt, G. M., Redding, S. J., Sturm, D. M., & Wolf, N. (2015). The Economics of Density: Evidence from the Berlin Wall. Econometrica, 83(6), 2127-2189. https://doi.org/10.3982/ECTA10876
