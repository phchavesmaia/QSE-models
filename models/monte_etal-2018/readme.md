# Monte, Redding and Rossi-Hansberg (2018)

This repository replicates Seidel and Wickerath (2020) implementation of the Monte, Redding, and Rossi-Hansberg (2018) model in Julia. It is **largely** based on the Ahfeldt and Seidel (2024) replication toolkit.

As emphasized by Ahfeldt and Seidel (2024) in their original implementation, the descriptive analyses and counterfactual exercises are for didactic purposes and unrelated to the research papers. Thus, this Julia implementation replicates the toolkit but does not reproduce results from either study.

Another important aspect to keep in mind is that my coding strategy prioritizes readability instead of efficiency. This is reflected by my active choice of element-wise operations instead of matricial operations, leading to higher computing times (when compared to the Matlab's implementation) despite Julia's efficient handling. In my honest opinion, the trade-off is very much worth it since the code remains rather fast (the slowest operation will less than 3 minutes to run) while being much more legible.

Thus, the main differences between my code and Ahfeldt and Seidel (2024) version are:

1. I consistently use element-wise operations, summing over a specific matrix dimension even if opting for a matricial multiplication implied higher efficiency;

2. In `solveProductTrade.jl` I define the tradable goods price index as:
  ```math
  P_n = \frac{\sigma}{\sigma - 1}\left( \frac{L_n^{1 - (1 - \sigma)\nu}}{\sigma f \pi_{nn}} \right)^{\frac{1}{1 - \sigma}}\frac{w_n d_{n n}}{A_{n}} 
  ```
  That is, as equation 11 of Seidel and Wickerath (2020). The author's of the original implementation define it, instead, as:
  ```math
  P_n = \frac{\sigma}{\sigma - 1}\left( \frac{L_n}{\sigma f \pi_{n n}} \right)^{\frac{1}{1 - \sigma}} \frac{w_n}{A_{n}}
  ```
  Which is equation 8 of Monte, Redding, and Rossi-Hansberg (2018) under the usual assumption that $`d_{nn}=1 \forall n\in\mathcal{N}`$.

3. While updating the value of $`\hat{P}_n`$, I allow for the (highly unlikely) possibility that $`\hat{d}_{nn}\neq 1`$, differently from the original authors, as stated in their Codebook. 
## References

G. M. Ahlfeldt, Seidel, T. (2024). Toolkit for quantitative spatial models. https://github.com/Ahlfeldt/MRRH2018-toolkit.

Monte, F., Redding, S. J., & Rossi-Hansberg, E. (2018). Commuting, migration, and local employment elasticities. American Economic Review, 108(12), 3855-3890. https://doi.org/10.1257/aer.20151507

Seidel, T., & Wickerath, J. (2020). Rush hours and urbanization. Regional Science and Urban Economics, 85, 103580. https://doi.org/10.1016/j.regsciurbeco.2020.103580v
