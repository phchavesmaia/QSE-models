# Quantitative Spatial Economics Models in Julia

These are my attempts at running QSE models in Julia. The goal is to discuss the steps and decisions made by the authors while replicating their findings.

Importantly, I value code readability more than computational efficiency, particularly for didactic purposes. Thus, I will consistently avoid turning the equilibrium equations into matrix multiplications, opting for element-wise multiplications instead. Nevertheless, Julia is very efficient when dealing with such operations, hence the computing times remain quite low across the repository.

It follows the roadmap of this repository:
- [X] [Redding and Rossi-Hansberg (2017)](https://github.com/phchavesmaia/QSE-models/tree/main/models/redding_rossihansberg-2017)
- [X] [Monte, Redding and Rossi-Hansberg (2018)](https://github.com/phchavesmaia/QSE-models/tree/main/models/monte_etal-2018)
- [ ] [Ahfeldt et al. (2015)](https://github.com/phchavesmaia/QSE-models/tree/main/models/ahfeldt_etal-2015)
- [ ] Hörcher and Graham (2024)

# Dependencies
If you are in a Windows machine, I highly recommend that you `instantiate` the repository as it is (Julia LTS). If you are in a Unix machine, it may be helpful to update `Manifest.toml` by accessing the Package Manager mode (`]`) and typing `resolve` followed by `update`.