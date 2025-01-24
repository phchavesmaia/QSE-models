This is my attempt to replicate Redding and Rossi-Hansberg (2017) implementation of Helpman (1998) in Julia. Notably, the biggest changes I have made in the code logic were:

1. Create the distance matrix using a graph (network) structure instead of the matricial distances
2. Rewrote a few of arguments that didn't made a lot of sense to me

I found that, apparently, Matlab reshapes matrices in a different way from Julia. Thus, straight up translating the code would imply running the model between a country at the North and the other at the South of the map, differently from the West-East division of the original paper. As of now, the code estimates $w_i$, $L_i$, and $\pi_{ni}$ in this inverted setting.

I am currently in the process of:
- [ ] Finish the "upside down" model
- [ ] "Rotate" the model so that it makes perfect sense in Julia logic
- [ ] Stress test the model to asses whether some normalizations (that I did not understand the reason behind) are actually necessary.

