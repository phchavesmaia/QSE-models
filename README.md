# Quantitative Spatial Economics Models in Julia

These are my attempts at running QSE models in Julia. This repository is mostly aspirational, since I added quite a lot of papers that employ this methodology and is highly unlikely that I (alone) will ever replicate them all. The goal is to discuss the steps and decisions made by the authors while replicating their findings. 

Importantly, since my primary goal is readability rather than efficiency, I will minimize the use of matrix multiplications, as they are often harder to read than simple element-wise (broadcasting) operations. In other words, I aim to prioritize Julia's ease for coding legible mathematical equations over its raw computational efficiency. 

It follows the "roadmap" of this repository:
- [X] [Redding and Rossi-Hansberg (2017)](https://github.com/phchavesmaia/QSE-models/tree/main/models/redding_rossihansberg-2017)
- [X] [Monte, Redding and Rossi-Hansberg (2018)](https://github.com/phchavesmaia/QSE-models/tree/main/models/monte_etal-2018)
- [ ] [Ahfeldt et al. (2015)](https://github.com/phchavesmaia/QSE-models/tree/main/models/ahfeldt_etal-2015)
- [ ] Hörcher and Graham (2024)
