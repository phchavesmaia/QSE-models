# Quantitative Spatial Economics Models in Julia

These are my attempts at running QSE models in Julia. This repository is mostly aspirational, since I added quite a lot of papers that employ this methodology and is highly unlikely that I (alone) will ever replicate them all. The goal is to discuss the steps and decisions made by the authors while replicating their findings.

Importantly, I value code readability more than computational efficiency, particularly for didactic purposes. Thus, I will consistently avoid turning the equilibrium equations into matrix multiplications, opting for element-wise multiplications instead. Nevertheless, Julia is very efficient when dealing with such operations, hence the computing times remain quite low across the repositories.

It follows the roadmap of this repository:
- [X] [Redding and Rossi-Hansberg (2017)](https://github.com/phchavesmaia/QSE-models/tree/main/models/redding_rossihansberg-2017)
- [X] [Monte, Redding and Rossi-Hansberg (2018)](https://github.com/phchavesmaia/QSE-models/tree/main/models/monte_etal-2018)
- [ ] [Ahfeldt et al. (2015)](https://github.com/phchavesmaia/QSE-models/tree/main/models/ahfeldt_etal-2015)
- [ ] Hörcher and Graham (2024)
