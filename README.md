# ISTA, FISTA and IRLS 

A repository of algorithms for sparse inversion for GEOPH-531 (Geophysical Inverse Problems). The idea is to obtain the sparse solution to the
problem 

$$A x - y \approx 0$$ 

which translates into minimizing the following cost function 

$$J = \|A x - y \|_2^2 + \lambda \|x \|_1$$

The notebook shows how to minimize J via ISTA, FISTA and IRLS, see the [PDF Document](https://github.com/msacchi/ISTA-FISTA-and-IRLS/blob/main/notes_ista.pdf) for a more detailed description of the problem. 

