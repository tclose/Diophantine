"""
Diophantine is a python package for finding small solutions of systems of
diophantine equations (see https://en.wikipedia.org/wiki/Diophantine_equation).
It is based on  PHP code by Keith Matthews (see www.number-theory.org) that
implements the algorithm described in the included 'algorithm.pdf' (see
http://www.numbertheory.org/lll.html for a list of associated publications).

There are two branches of this code in the GitHub repository
(see https://github.com/tclose/Diophantine.git), 'master', which uses the
sympy library and therefore uses arbitrarily long integer representations, and
'numpy', which uses the numpy library, which is faster but can suffer from
integer overflow errors despite using int64 representations.

Diophantine is released under the MIT Licence (see Licence for details)

Author: Thomas G. Close (tom.g.close@gmail.com)
"""
# The MIT License (MIT)
#
# Copyright (c) 2015 Thomas G. Close
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
from copy import deepcopy
from math import copysign, sqrt, log10, floor
from fractions import gcd
from sympy import Matrix, zeros, ones, eye
from itertools import chain, product
global print_count
print_count = 0
verbose_solve = False
verbose_hnf = False
verbose_chol = False


class NoSolutionException(Exception):
    pass


# Sign of a variable, which isn't included in math for some reason
def sign(x):
    return copysign(1, x) if x else 0


def nonzero(m):
    return [(i, j) for i, j in product(xrange(m.shape[0]), xrange(m.shape[1]))
            if m[i, j] != 0]


def solve(A, b):
    """
    Finds small solutions to systems of diophantine equations, A x = b, where A
    is a M x N matrix of coefficents, b is a M x 1 vector and x is the
    N x 1 solution vector, e.g.

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

    The returned solution vector will tend to be one with the smallest norms.
    If multiple solutions with the same norm are found they will all be
    returned. If there are no solutions the empty list will be returned.
    """
    A = Matrix(A)
    b = Matrix(b)
    if b.shape != (A.shape[0], 1):
        raise Exception("Length of b vector ({}) does not match number of rows"
                        " in A matrix ({})".format(b.shape[0], A.shape[0]))
    if verbose_solve:
        Ab = A.row_join(b)
        print 'Ab: "' + ' '.join(str(v) for v in Ab.T.vec()) + '"'
    G = zeros(A.shape[1] + 1, A.shape[0] + 1)
    G[:-1, :-1] = A.T
    G[-1, :-1] = b.reshape(1, b.shape[0])
    G[-1, -1] = 1
    # A is m x n, b is m x 1, solving AX=b, X is n x 1+
    # Ab is the (n+1) x m transposed augmented matrix. G=[A^t|0] [b^t]1]
    if verbose_solve:
        print "G:"
        printnp(G)
    hnf, P, rank = lllhermite(G)
    if verbose_solve:
        print "HNF(G):"
        printnp(hnf)
        print "P:"
        printnp(P)
        print "Rank: {}".format(rank)
    r = rank - 1  # For convenience
    if not any(chain(hnf[:r, -1], hnf[r, :-1])) and hnf[r, -1] == 1:
        nullity = hnf.shape[0] - rank
        if nullity:
            basis = P[rank:, :-1].col_join(-P[r, :-1])
            if verbose_solve:
                print "Basis:\n"
                printnp(basis)
            solutions = get_solutions(basis)
        else:
            raise NotImplementedError("Ax=B has unique solution in integers")
    else:
        solutions = []
    return solutions


def lllhermite(G, m1=1, n1=1):
    """
    Input: integer mxn matrix A, nonzero, at least two rows+
    Output: small unimodular matrix B and HNF(A), such that BA=HNF(A)+
    The Havas, Majewski, Matthews LLL method is used+
    We usually take alpha=m1/n1, with (m1,n1)=(1,1) to get best results+
    """
    m = G.shape[0]
    n = G.shape[1]
    A, B, L, D = initialise_working_matrices(G)
    if first_nonzero_is_negative(A):
        B[m, m] = -1
        A[m, :] *= -1
    k = 1
    while k < m:
        if verbose_hnf:
            print "k={k}, m={m}".format(k=k, m=m)
        col1, col2 = reduce_matrix(A, B, L, k, k - 1, D)
        if verbose_hnf:
            print "col1={col1}, col2={col2}".format(col1=col1, col2=col2)
        if verbose_hnf:
            print_all(A, B, L, D)
        u = n1 * (int(D[k - 1]) * int(D[k + 1]) +
                  int(L[k, k - 1]) * int(L[k, k - 1]))
        v = m1 * int(D[k]) * int(D[k])
        if verbose_hnf:
            print "u={u}, v={v}".format(u=u, v=v)
        if col1 <= min(col2, n - 1) or (col1 == n and col2 == n and u < v):
            swap_rows(k, A, B, L, D)
            if verbose_hnf:
                print_all(A, B, L, D)
            if k > 1:
                k = k - 1
                if verbose_hnf:
                    print "col1 <= minim && k > 1"
            else:
                if verbose_hnf:
                    print "col1 <= minim"
        else:
            for i in reversed(xrange(k - 1)):
                reduce_matrix(A, B, L, k, i, D)
            k = k + 1
            if verbose_hnf:
                print "col1 > minim"
    try:
        rank = A.shape[0] - next(i for i in xrange(A.shape[0])
                                 if nonzero(A[i, :]))
    except StopIteration:
        assert False, "A matrix contains only zeros"
    hnf = A[::-1, :]
    unimodular_matrix = B[::-1, :]
    return hnf, unimodular_matrix, rank


def initialise_working_matrices(G):
    """  G is a nonzero matrix with at least two rows.  """
    B = eye(G.shape[0])
    # Lower triang matrix
    L = zeros(G.shape[0], G.shape[0])
    D = ones(G.shape[0] + 1, 1)
    A = Matrix(G)
    return A, B, L, D


def first_nonzero_is_negative(A):
    """
    returns 0 if the first nonzero column j of A contains more than one nonzero
    entry, or contains only one nonzero entry and which is positive+ returns 1
    if the first nonzero column j of A contains only one nonzero entry, which
    is negative+ This assumes A is a nonzero matrix with at least two rows+
    """
    nonzero_columns = zip(*nonzero(A))[1]  # Should always be nonzero
    # Get the first nonzero column
    nonzero_col = A[:, min(nonzero_columns)]
    # Get the nonzero elements
    nonzero_elems = [e for e in nonzero_col if e != 0]
    # If there is only one and it is negative return 1 else 0
    return len(nonzero_elems) == 1 and nonzero_elems[0] < 0


def reduce_matrix(A, B, L, k, i, D):
    nonzero_i_elems = zip(*nonzero(A[i, :]))
    if len(nonzero_i_elems):
        col1 = nonzero_i_elems[1][0]
        if A[i, col1] < 0:
            minus(i, L)
            A[i, :] *= -1
            B[i, :] *= -1
    else:
        col1 = A.shape[1]
    nonzero_k_elems = zip(*nonzero(A[k, :]))
    if len(nonzero_k_elems):
        col2 = nonzero_k_elems[1][0]
    else:
        col2 = A.shape[1]
    if col1 < A.shape[1]:
        q = A[k, col1] // A[i, col1]
    else:
        t = abs(L[k, i])
        t = 2 * t
        if t > D[i + 1]:
            q = lnearint(L[k, i], D[i + 1])
        else:
            q = 0
    if q != 0:
        A[k, :] -= q * A[i, :]
        B[k, :] -= q * B[i, :]
        L[k, i] -= q * D[i + 1]
        L[k, :i] -= q * L[i, :i]
    return col1, col2


def minus(j, L):
    L[j, :] = -L[j, :]
    L[:, j] = -L[:, j]


def swap_rows(k, A, B, L, D):
    # To avoid the interpretation of -1 as the last index of the matrix create
    # a reverse stop that ends past the negative of the length of the matrix
    reverse_stop = k - 2 if k > 1 else -(A.shape[0] + 1)
    # Swap rows of the matrices
    A[(k - 1):(k + 1), :] = A[k:reverse_stop:-1, :]
    B[(k - 1):(k + 1), :] = B[k:reverse_stop:-1, :]
    L[(k - 1):(k + 1), :(k - 1)] = L[k:reverse_stop:-1, :(k - 1)]
    if verbose_hnf:
        print_all(A, B, L, D)
    t = (L[(k + 1):, k - 1] * D[k + 1] / D[k] -
         L[(k + 1):, k] * L[k, k - 1] / D[k])
    L[(k + 1):, k - 1] = (L[(k + 1):, k - 1] * L[k, k - 1] +
                          L[(k + 1):, k] * D[k - 1]) / D[k]
    L[(k + 1):, k] = t
    if verbose_hnf:
        print_all(A, B, L, D)
    t = int(D[k - 1]) * int(D[k + 1]) + int(L[k, k - 1]) * int(L[k, k - 1])
    D[k] = t / D[k]


def get_solutions(A):
    m = A.shape[0]
    n = A.shape[1]
    G = gram(A)
    if verbose_solve:
        print "G:"
        printnp(G)
    N, D = cholesky(G)
    Qn, Qd = N, D
    if verbose_solve:
        print "Qn:"
        printnp(Qn)
        print "Qd:"
        printnp(Qd)
    m -= 1
    Nn = Qn[:m, m]
    Nd = Qd[:m, m]
    Cn = 0
    Cd = 1
    if verbose_solve:
        print "N:"
        printnp(Nn)
        print "D:"
        printnp(Nd)
    for i in xrange(m):
        num, den = multr(Nn[i], Nd[i], Nn[i], Nd[i])
        num, den = multr(num, den, Qn[i, i], Qd[i, i])
        Cn, Cd = addr(Cn, Cd, num, den)
        if verbose_solve:
            print "i: {}, Cnum: {}, Cden: {}".format(i + 1, Cn, Cd)
    i = m - 1
    # List to hold working variables
    x = zeros(m, 1)
    UB = zeros(m, 1)
    Tn = zeros(m, 1)
    Td = zeros(m, 1)
    Un = zeros(m, 1)
    Ud = zeros(m, 1)
    Tn[i] = Cn
    Td[i] = Cd
    Un[i] = 0
    Ud[i] = 1
    solutions = []  # List to hold multipliers
    while True:
        # Calculate UB
        Zn, Zd = ratior(Tn[i], Td[i], Qn[i, i], Qd[i, i])
        num, den = subr(Nn[i], Nd[i], Un[i], Ud[i])
        if verbose_solve:
            print "Tn:"
            printnp(Tn[i:])
            print "Td:"
            printnp(Td[i:])
            print "Zn: {}, Zd: {}, num: {}, den: {}".format(
                Zn, Zd, num, den)
        UB[i] = introot(Zn, Zd, num, den)
        # Calculate x
        num, den = subr(Un[i], Ud[i], Nn[i], Nd[i])
        x[i] = -introot(Zn, Zd, num, den) - 1
        while True:
            if verbose_solve:
                print "i: {}, UB: {}".format(i, UB[i])
                print "x:"
                printnp(x[i:])
            x[i] += 1
            if x[i] <= UB[i]:
                if i == 0:
                    if verbose_solve:
                        print "x:"
                        printnp(x[i:])
                    lcv = lcasvector(A[:-1, :], x)
                    if verbose_solve:
                        print "lcv:"
                        printnp(lcv)
                    solution = A[m, :n] - lcv.reshape(1, lcv.shape[0])
                    if verbose_solve:
                        print "solution:"
                        printnp(solution)
                    solutions.append(solution.T)
                else:
                    # now update U
                    Un[i - 1], Ud[i - 1] = 0, 1
                    for j in xrange(i, m):
                        # Loops from back of xs
                        num, den = multr(Qn[i - 1, j], Qd[i - 1, j], x[j], 1)
                        Un[i - 1], Ud[i - 1] = addr(Un[i - 1], Ud[i - 1], num,
                                                    den)
                        if verbose_solve:
                            print ("i: {}, j: {}, Un: {}, Ud: {}, num: {}, "
                                   "den: {}".format(i, j, Un[i - 1], Ud[i - 1],
                                                    num, den))
                    # now update T
                    num, den = addr(x[i], 1, Un[i], Ud[i])
                    num, den = subr(num, den, Nn[i], Nd[i])
                    num, den = multr(num, den, num, den)
                    num, den = multr(Qn[i, i], Qd[i, i], num, den)
                    Tn[i - 1], Td[i - 1] = subr(Tn[i], Td[i], num, den)
                    i = i - 1
                    break
            else:
                i = i + 1
                if i == m:
                    return solutions


def cholesky(A):
    """
    # A is positive definite mxm
    """
    assert A.shape[0] == A.shape[1]
    # assert all(A.eigenvals() > 0)
    m = A.shape[0]
    N = deepcopy(A)
    D = ones(*A.shape)
    for i in xrange(m - 1):
        for j in xrange(i + 1, m):
            N[j, i] = N[i, j]
            D[j, i] = D[i, j]
            n, d = ratior(N[i, j], D[i, j], N[i, i], D[i, i])
            N[i, j], D[i, j] = n, d
            if verbose_chol:
                print "i={}, j={}".format(i + 1, j + 1)
                print "N:"
                printnp(N)
                print "D:"
                printnp(D)
        for k in xrange(i + 1, m):
            for l in xrange(k, m):
                n, d = multr(N[k, i], D[k, i], N[i, l], D[i, l])
                N[k, l], D[k, l] = subr(N[k, l], D[k, l], n, d)
                if verbose_chol:
                    print "k={}, l={}".format(k + 1, l + 1)
                    print "N:"
                    printnp(N)
                    print "D:"
                    printnp(D)
    return N, D


def gram(A):
    """
    Need to check for row and column operations
    """
    m = A.shape[0]
    B = zeros(m, m)
    for i in xrange(m):
        for j in xrange(m):
            B[i, j] = A[i, :].dot(A[j, :])  # dotproduct(A[i], A[j], n)
    return Matrix(B)


def introot(a, b, c, d):
    """
    With Z=a/b, U=c/d, returns [sqrt(a/b)+c/d]. First ANSWER =
    [sqrt(Z)] + [U]. One then tests if Z < ([sqrt(Z)] + 1 -U)^2. If
    this does not hold, ANSWER += 1+ For use in fincke_pohst()+
    """
    y = c // d
    if a == 0:
        return y
    x = a // b
    assert x >= 0
    x_sqrt = int(floor(sqrt(x)))
    answer = x_sqrt + y
    num, den = subr(c, d, y, 1)
    num, den = subr(1, 1, num, den)
    num, den = addr(x_sqrt, 1, num, den)
    num, den = multr(num, den, num, den)
    t = comparer(num, den, a, b)
    if t <= 0:
        answer = answer + 1
    int_answer = int(answer)
    assert int_answer == answer
    return int_answer


def egcd(p, q):
    if q == 0:
        if p != 0:
            s = sign(p)
            if s == 1:
                k1 = 1
            else:
                k1 = -1
            return abs(p), k1, 0
        else:
            return 0, 0, 0
    a = p
    b = abs(q)
    c = a % b
    s = sign(q)
    if c == 0:
        if s == 1:
            k2 = 1
        else:
            k2 = -1
        return b, 0, k2
    l1 = 1
    k1 = 0
    l2 = 0
    k2 = 1
    while c != 0:
        q = a // b
        a = b
        b = c
        c = a % b
        h1 = l1 - q * k1
        h2 = l2 - q * k2
        l1 = k1
        l2 = k2
        k1 = h1
        k2 = h2
    if s == -1:
        k2 = 0 - k2
    return b, k1, k2


def lnearint(a, b):
    """
    left nearest integer
    returns y+1/2 if a/b=y+1/2, y integral+
    """
    y = a // b
    if b < 0:
        a = -a
        b = -b
    x = b * y
    z = a - x
    z = 2 * z
    if z > b:
        y = y + 1
    return y


def ratior(a, b, c, d):
    """ returns (a/b)/(c/d)"""
    r = a * d
    s = b * c
    g = abs(gcd(r, s))
    if s < 0:
        g = -g
    return r / g, s / g


def multr(a, b, c, d):
    # returns (a/b)(c/d)
    r = a * c
    s = b * d
    g = abs(gcd(r, s))
    return r / g, s / g


def subr(a, b, c, d):
    t = a * d - b * c
    u = b * d
    g = abs(gcd(t, u))
    return t / g, u / g


def addr(a, b, c, d):
    t = a * d + b * c
    u = b * d
    g = abs(gcd(t, u))
    return t / g, u / g


def comparer(a, b, c, d):
    """Assumes b>0 and d>0.  Returns -1, 0 or 1 according as a/b <,=,> c/d+ """
    assert b > 0 and d > 0
    return sign(a * d - b * c)


def lcasvector(A, x):
    """lcv[j]=X[1]A[1, j]=...+X[m]A[m, j], 1 <= j <= n+"""
    # global lcv
#     printnp(x)
#     printnp(A)
    n = A.shape[1]
    lcv = zeros(n, 1)
    for j in xrange(n):
        lcv[j] = x.dot(A[:, j])
    return lcv


def print_all(A, B, L, D):
    global print_count
    print "------ print {} -----".format(print_count)
    print 'A: '
    printnp(A)
    print 'B: '
    printnp(B)
    print 'L: '
    printnp(L)
    print 'D: '
    printnp(D)
    print_count += 1


def num_chars(d):
    if d == 0:
        num = 2
    else:
        num = int(floor(log10(abs(d)))) + int(sign(d) < 0) + 2
    return num


def printnp(m):
    """
    Prints a sympy matrix in the same format as a numpy array
    """
    m = Matrix(m)
    if m.shape[1] == 1:
        m = m.reshape(1, m.shape[0])
    elem_len = max(num_chars(d) for d in m.vec())
    if m.shape[0] > 1:
        outstr = '[['
    else:
        outstr = '['
    for i in xrange(m.shape[0]):
        if i:
            outstr += ' ['
        char_count = 0
        for j, elem in enumerate(m[i, :]):
            char_count += elem_len
            if char_count > 77:
                outstr += '\n '
                char_count = elem_len + 2
            # Add spaces
            outstr += ' ' * (elem_len - num_chars(elem) +
                             int(j != 0)) + str(elem)

        if i < m.shape[0] - 1:
            outstr += ']\n'
    if m.shape[0] > 1:
        outstr += ']]'
    else:
        outstr += ']'
    print outstr
