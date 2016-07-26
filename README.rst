Diophantine
===========

Diophantine (http://github.com/tclose/Diophantine) is a Python package for
finding small (integer) solutions of systems of diophantine equations
(see http://en.wikipedia.org/wiki/Diophantine_equation). It is based on  PHP
code by Keith Matthews (see www.number-theory.org) that implements the
algorithm described in https://github.com/tclose/Diophantine/blob/master/algorithm.pdf
(see http://www.numbertheory.org/lll.html for a list of associated publications),
which uses the LLL algorithm to calculate the Hermite-normal-form described in
the paper:

Extended gcd and Hermite normal form algorithms via lattice basis reduction,
G. Havas, B.S. Majewski, K.R. Matthews, Experimental Mathematics, Vol 7 (1998) 125-136

(please cite this paper if you use this code in a scientific publication)

There are two branches of this code in the GitHub repository 
(see https://github.com/tclose/Diophantine.git), 'master', which uses the
sympy library and therefore uses arbitrarily long integer representations, and 
'numpy', which uses the numpy library, which is faster but can suffer from
integer overflow errors despite using int64 representations

To find small solutions to a system of diophantine equations, A x = b, where A
is a M x N matrix of coefficents, b is a M x 1 vector and x is the
N x 1 vector, use the 'solve' method in the module, e.g.

    >>> from sympy import Matrix
    >>> from diophantine import solve
    >>> A = Matrix([[1, 0, 0, 2], [0, 2, 3, 5], [2, 0, 3, 1], [-6, -1, 0, 2],
                    [0, 1, 1, 1], [-1, 2, 0,1], [-1, -2, 1, 0]]).T
    >>> b = Matrix([1, 1, 1, 1])
    >>> solve(A, b)
    [Matrix([
    [-1],
    [ 1],
    [ 0],
    [ 0],
    [-1],
    [-1],
    [-1]])]

The returned solution vector will tend to be one with the smallest norms. If multiple solutions with the same norm are found they will all be returned. If there are no solutions the empty list will be returned.

Diophantine is released under the MIT Licence (see Licence for details)

Author: Thomas G. Close (tom.g.close@gmail.com)

.. image:: https://travis-ci.org/tclose/Diophantine.svg?branch=master
   :target: https://travis-ci.org/tclose/Diophantine
   :alt: Unit Test Status
.. image:: https://coveralls.io/repos/github/tclose/Diophantine/badge.svg?branch=master
   :target: https://coveralls.io/github/tclose/Diophantine?branch=master
   :alt: Unit Test Coverage
