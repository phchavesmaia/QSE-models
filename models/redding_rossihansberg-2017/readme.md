# Redding and Rossi-Hansberg (2017) 

This is my attempt to replicate Redding and Rossi-Hansberg (2017) implementation of the Helpman (1998) model in Julia. Notably, the biggest changes I have made in the code logic were:

1. Create the distance matrix using a graph (network) structure instead of the matricial distances
2. Rewrote a couple of functions that didn't made a lot of sense to me, though they still lead to identical results

I found that, apparently, Matlab reshapes matrices in a different way from Julia. Thus, straight up translating the code would imply running the model for a North-South country split, differently from the West-East division of the original paper. In this implementation of the code, this incongruity between languages is fixed, i.e., the code replicates the results of the paper, although the matrices become slighly different.

## References

Helpman E. (1998). The size of regions. In Topics in Public Economics: Theoretical and Applied Analysis, ed. D Pines, E Sadka, I Zilcha, pp. 33–54. Cambridge, UK: Cambridge Univ. Press. (Working paper version: https://matthewturner.org/ec2410/readings/Helpman_unp_1995.pdf )

Redding, S. J., & Rossi-Hansberg, E. (2017). Quantitative spatial economics. Annual Review of Economics, 9(1), 21-58. https://doi.org/10.1146/annurev-economics-063016-103713

---

I owe quite a lot to [@SebKrantz](https://github.com/SebKrantz/Quantitative-Spatial-Economics/tree/main/QSE-ARE-2017) who also took upon himself to accomplish this challange.
