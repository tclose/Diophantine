Diophantine
======

Diophantine is a python package for finding small solutions of systems of
diophantine equations (see https://en.wikipedia.org/wiki/Diophantine_equation).
It is based on  PHP code by Keith Matthews (see www.number-theory.org) that
implements the algorithm described in the included 'algorithm.pdf' (see
http://www.numbertheory.org/lll.html for a list of associated publications),
which uses the LLL algorith to calculate the Hermite-normal-form described in
the paper:

Extended gcd and Hermite normal form algorithms via lattice basis reduction,
G. Havas, B.S. Majewski, K.R. Matthews, Experimental Mathematics, Vol 7 (1998) 125-136.

(please cite this paper if you use this code in a scientific publication)

There are two branches of this code in the GitHub repository 
(see https://github.com/tclose/Diophantine.git), 'master', which uses the
sympy library and therefore uses arbitrarily long integer representations, and 
'numpy', which uses the numpy library, which is faster but can suffer from
integer overflow errors despite using int64 representations.

Diophantine is released under the MIT Licence (see Licence for details)

Author: Thomas G. Close (tom.g.close@gmail.com)
