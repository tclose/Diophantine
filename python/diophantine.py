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

# Original written in PHP by Keith Matthews (webmaster@number-theory.org)
#
# Converted from PHP to Python by Thomas G. Close (tom.g.close@gmail.com)
import os
from copy import deepcopy
from fractions import gcd
import numpy
from nineml import units as un
from itertools import chain
global print_count
print_count = 0

numpy.seterr(over='raise')
# import warnings
# warnings.filterwarnings('error')


def solve(A, b):
    """
    Finds the minimal combination of reference dimensions to make the compound
    dimension
    """
    A = numpy.asarray(A, dtype=numpy.int64)
    b = numpy.asarray(b, dtype=numpy.int64)
    G = numpy.zeros((A.shape[1] + 1, A.shape[0] + 1), dtype=numpy.int64)
    G[:-1, :-1] = A.T
    G[-1, :-1] = b
    G[-1, -1] = 1
    # A is m x n, b is m x 1, solving AX=b, X is n x 1+
    # Ab is the (n+1) x m transposed augmented matrix. G=[A^t|0] [b^t]1]
    print "G:"
    print G
    hnf, P, rank = lllhermite(G)
    print "HNF(G):"
    print hnf
    print "P:"
    print P
    print "Rank: {}".format(rank)
    r = rank - 1  # For convenience
    if not any(chain(hnf[:r, -1], hnf[r, :-1])) and hnf[r, -1] == 1:
        nullity = hnf.shape[0] - rank
        if nullity:
            basis = numpy.concatenate((P[rank:, :-1],
                                       -P[r, :-1].reshape(1, -1)))
            print "Basis:\n"
            print basis
            x = shortest_distance_axb(basis)
        else:
            raise NotImplementedError("Ax=B has unique solution in integers")
    else:
        raise Exception("AX=B has no solution in integers")
    return x


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
        #print "k={k}, m={m}".format(k=k, m=m)
        col1, col2 = reduce_matrix(A, B, L, k, k - 1, D)
        #print "col1={col1}, col2={col2}".format(col1=col1, col2=col2)
        #print_all(A, B, L, D)
        u = n1 * (int(D[k - 1]) * int(D[k + 1]) +
                  int(L[k, k - 1]) * int(L[k, k - 1]))
        v = m1 * int(D[k]) * int(D[k])
        #print "u={u}, v={v}".format(u=u, v=v)
        if col1 <= min(col2, n - 1) or (col1 == n and col2 == n and u < v):
            swap_rows(k, A, B, L, D)
            #print_all(A, B, L, D)
            if k > 1:
                k = k - 1
                #print "col1 <= minim && k > 1"
            else:
                #print "col1 <= minim"
                pass
        else:
            for i in reversed(xrange(k - 1)):
                reduce_matrix(A, B, L, k, i, D)
            k = k + 1
            #print "col1 > minim"
    try:
        rank = len(A) - next(i for i, row in enumerate(A) if any(row != 0))
    except StopIteration:
        assert False, "A matrix contains only zeros"
    hnf = A[::-1, :]
    unimodular_matrix = B[::-1, :]
    return hnf, unimodular_matrix, rank


def initialise_working_matrices(G):
    """  G is a nonzero matrix with at least two rows.  """
    B = numpy.eye(G.shape[0], dtype=numpy.int64)
    # Lower triang matrix
    L = numpy.zeros((G.shape[0], G.shape[0]), dtype=numpy.int64)
    D = numpy.ones(G.shape[0] + 1, dtype=numpy.int64)
    A = numpy.array(G, dtype=numpy.int64)
    return A, B, L, D


def first_nonzero_is_negative(A):
    """
    returns 0 if the first nonzero column j of A contains more than one nonzero
    entry, or contains only one nonzero entry and which is positive+ returns 1
    if the first nonzero column j of A contains only one nonzero entry, which
    is negative+ This assumes A is a nonzero matrix with at least two rows+
    """
    nonzero_columns = numpy.nonzero(
        numpy.sum(A, axis=0, dtype=numpy.int64) != 0)[0]
    assert len(nonzero_columns)
    # Get the first nonzero column
    nonzero_col = A[:, numpy.min(nonzero_columns)]
    # Get the nonzero elements
    nonzero_elems = numpy.nonzero(nonzero_col)[0]
    # If there is only one and it is negative return 1 else 0
    return len(nonzero_elems) == 1 and nonzero_elems[0] < 0


def reduce_matrix(A, B, L, k, i, D):
    nonzero_i_elems = numpy.nonzero(A[i])[0]
    if len(nonzero_i_elems):
        col1 = nonzero_i_elems[0]
        if A[i, col1] < 0:
            minus(i, L)
            A[i, :] *= -1.0
            B[i, :] *= -1.0
    else:
        col1 = A.shape[1]
    nonzero_k_elems = numpy.nonzero(A[k])[0]
    if len(nonzero_k_elems):
        col2 = nonzero_k_elems[0]
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
    A[(k - 1, k), :] = A[(k, k - 1), :]
    B[(k - 1, k), :] = B[(k, k - 1), :]
    L[(k - 1, k), :(k - 1)] = L[(k, k - 1), :(k - 1)]
#     #print_all(A, B, L, D)
    t = (L[(k + 1):, k - 1] * D[k + 1] / D[k] -
         L[(k + 1):, k] * L[k, k - 1] / D[k])
    L[(k + 1):, k - 1] = (L[(k + 1):, k - 1] * L[k, k - 1] +
                          L[(k + 1):, k] * D[k - 1]) / D[k]
    L[(k + 1):, k] = t
#     #print_all(A, B, L, D)
    t = int(D[k - 1]) * int(D[k + 1]) + int(L[k, k - 1]) * int(L[k, k - 1])
    D[k] = t / D[k]


def shortest_distance_axb(A):
    m = A.shape[0]
    n = A.shape[1]
    G = gram(A)
    print "G:"
    print G
    N, D = cholesky(G)
    Qn, Qd = N, D
    print "Qn:"
    print Qn
    print "Qd:"
    print Qd
    m -= 1
    Nn = Qn[:m, m]
    Nd = Qd[:m, m]
    Cn = 0
    Cd = 1
    print "N:"
    print Nn
    print "D:"
    print Nd
    for i in xrange(m):
        num, den = multr(Nn[i], Nd[i], Nn[i], Nd[i])
        num, den = multr(num, den, Qn[i][i], Qd[i][i])
        Cn, Cd = addr(Cn, Cd, num, den)
    i = m - 1
    Tn = Cn
    Td = Cd
    Un = 0
    Ud = 1
    solutions = []  # List to hold multipliers
    x = numpy.empty(m, dtype=numpy.int64)  # List to hold pas values of x
    while 1:
        # Calculate UB
        Zn, Zd = ratior(Tn, Td, Qn[i, i], Qd[i, i])
        num, den = subr(Nn[i], Nd[i], Un, Ud)
        UB = introot(Zn, Zd, num, den)
        # Calculate x
        num, den = subr(Un, Ud, Nn[i], Nd[i])
        x[i] = -introot(Zn, Zd, num, den) - 1
        print "x:"
        print x[i:]
        while True:
            x[i] += 1
            if x[i] <= UB:
                if i == 0:
                    print "x:"
                    print x[i:]
                    lcv = lcasvector(A[:-1, :], x)
                    print "lcv:"
                    print lcv
                    solution = A[m, :n] - lcv
                    print "solution:"
                    print solution
                    solutions.append(solution)
                    continue
                else:
                    # now update U
                    prev_Un = Un
                    prev_Ud = Ud
                    Un, Ud = 0, 1
                    for j in xrange(i + 1, m):
                        # Loops from back of xs
                        num, den = multr(Qn[i, j], Qd[i, j], x[j], 1)
                        Un, Ud = addr(Un, Ud, num, den)
                    # now update T
                    num, den = addr(x[i], 1, prev_Un, prev_Ud)
                    num, den = subr(num, den, Nn[i], Nd[i])
                    num, den = multr(num, den, num, den)
                    num, den = multr(Qn[i][i], Qd[i][i], num, den)
                    Tn, Td = subr(Tn, Td, num, den)
                    i = i - 1
                    break
            else:
                i = i + 1
                if i == m:
                    return solutions
                continue


def cholesky(A):
    """
    # A is positive definite mxm
    """
    assert A.ndim == 2 and A.shape[0] == A.shape[1]
    assert numpy.all(numpy.linalg.eigvals(A) > 0)
    m = A.shape[0]
    N = deepcopy(A)
    D = numpy.ones(A.shape, dtype=numpy.int64)
    for i in xrange(m - 1):
        for j in xrange(i + 1, m):
            N[j][i] = N[i][j]
            D[j][i] = D[i][j]
            n, d = ratior(N[i][j], D[i][j], N[i][i], D[i][i])
            N[i][j], D[i][j] = n, d
            print "i={}, j={}".format(i + 1, j + 1)
            print "N:"
            print N
            print "D:"
            print D
        for k in xrange(i + 1, m):
            for l in xrange(k, m):
                n, d = multr(N[k][i], D[k][i], N[i][l], D[i][l])
                N[k][l], D[k][l] = subr(N[k][l], D[k][l], n, d)
                print "k={}, l={}".format(k + 1, l + 1)
                print "N:"
                print N
                print "D:"
                print D
    return N, D


def gram(A):
    """
    Need to check for row and column operations
    """
    m = A.shape[0]
    B = numpy.empty((m, m), dtype=numpy.int64)
    for i in xrange(m):
        for j in xrange(m):
            B[i][j] = A[i].dot(A[j])  # dotproduct(A[i], A[j], n)
    return numpy.array(B, dtype=numpy.int64)


def introot(a, b, c, d):
    """
    With Z=a/b, U=c/d, returns [numpy.sqrt(a/b)+c/d]. First ANSWER =
    [numpy.sqrt(Z)] + [U]. One then tests if Z < ([numpy.sqrt(Z)] + 1 -U)^2. If
    this does not hold, ANSWER += 1+ For use in fincke_pohst()+
    """
    y = c // d
    if a == 0:
        return y
    x = a / b
    assert x >= 0
    x = numpy.sqrt(x)
    answer = x + y
    n, d = subr(c, d, y, 1)
    n, d = subr(1, 1, n, d)
    n, d = addr(x, 1, n, d)
    n, d = multr(n, d, n, d)
    t = comparer(n, d, a, b)
    if t <= 0:
        answer = answer + 1
    int_answer = int(answer)
    assert int_answer == answer
    return int_answer


def egcd(p, q):
    if q == 0:
        if p != 0:
            s = numpy.sign(p)
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
    s = numpy.sign(q)
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
    return numpy.sign(a * d - b * c)


def lcasvector(A, x):
    """lcv[j]=X[1]A[1][j]=...+X[m]A[m][j], 1 <= j <= n+"""
    # global lcv
#     print x
#     print A
    n = A.shape[1]
    lcv = numpy.empty(n, dtype=numpy.int64)
    for j in xrange(n):
        lcv[j] = x.dot(A[:, j])
    return lcv


if __name__ == '__main__':
    reference_dims = [un.current, un.time, un.voltage,
                      un.specificCapacitance,
                      un.conductanceDensity, un.luminous_intensity,
                      un.temperature, un.substance]
    compound = (un.voltage * un.temperature) / un.specificCapacitance
#     diophantine(compound, reference_dims)


arrays = [
    numpy.array([
        [-3, -2, -4, -3, -1, 0, -3, 0, 1, 3],
        [3, -4, 3, -1, 3, -2, -4, -2, -1, 0],
        [2, 1, 0, -2, -4, 3, 3, -4, 0, 0],
        [4, 4, 3, -4, 2, 4, 1, 0, -3, -2],
        [1, 2, 2, 1, -2, 0, 2, 0, -3, -1],
        [4, 0, -2, -1, 0, 4, 4, 2, 0, 0],
        [-4, 1, -4, 4, -4, 0, -2, 3, 4, 4]], dtype=numpy.int64),
    numpy.array([
        [-1, 0, 0, -3, -3, -3, 4, 0, 1, 4],
        [0, -2, -2, 4, 2, -4, 0, -3, -4, 2],
        [-2, 3, 1, -4, 2, -1, 1, -4, 0, 1],
        [4, -3, -2, 2, -1, 1, -4, -2, 4, 1],
        [-3, -2, -1, -3, 0, -4, 1, -3, 3, 1],
        [-4, -1, 0, -3, 0, 0, 3, 3, -4, 0],
        [1, 1, -1, -3, 2, 2, -3, 3, 2, 2]], dtype=numpy.int64),
    numpy.array([
        [-1, -2, 3, 0, 4, 0, -4, -3, 4, -2],
        [-3, 2, -2, 0, -4, 3, 3, 2, 0, -4],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, -1, 0, 0, -1, -1, 3, -3, 4, 2],
        [3, -1, 2, 0, -4, -3, -1, -1, 2, 3],
        [-4, 0, -1, 0, 4, 0, 1, -4, -2, 0],
        [4, 2, 3, 0, 0, -2, 2, -2, -4, 1]], dtype=numpy.int64),
    numpy.array([
        [-3, 3, 4, -1, 0, -4, -1, -4, 2, -2],
        [1, 2, 3, -1, -3, 3, -3, -2, 1, -2],
        [-4, 2, 2, -2, -3, -1, -2, -4, 0, 2],
        [0, -4, -3, -3, 1, 2, 0, -3, 1, -1],
        [-1, -1, 3, 1, 1, 4, -3, -3, 0, 2],
        [0, 1, -4, 1, -3, 0, -1, 0, 1, 0],
        [0, 0, -2, -2, 4, 0, 4, 1, 2, 0]], dtype=numpy.int64),
    numpy.array([
        [4, -3, 0, 0, 0, 3, 4, -4, 0, -3],
        [-4, 4, -3, 0, -3, -3, 2, 0, -1, -1],
        [0, 3, 4, 0, 2, -2, 2, 2, 0, 3],
        [-3, 1, 0, 0, 2, 0, 0, -3, 1, 1],
        [0, -4, -3, 0, 0, 1, -3, -1, 1, 0],
        [4, 3, 2, 0, 1, -1, 0, -2, 2, -2],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]], dtype=numpy.int64)]

xs = [
    numpy.array([3, 4, 1, 0, -3, 0, 1, 4, -2], dtype=numpy.int64),
    numpy.array([0, 4, 4, 3, 2, 4, 0, 1, 4], dtype=numpy.int64),
    numpy.array([0, -4, 0, 1, -4, 4, 4, -2, 3], dtype=numpy.int64),
    numpy.array([1, 2, -3, 3, -4, 1, -3, -3, -4], dtype=numpy.int64),
    numpy.array([0, -4, -1, -2, -4, 0, 4, 3, -4], dtype=numpy.int64)]


# for i, arr in enumerate(arrays):
#     print '$arrays[{}] = "{}";'.format(
#         i, " ".join([str(e) for e in arr.ravel()]))
# quit()

# for i, x in enumerate(xs):
#     print '$xs[{}] = "{}";'.format(
#         i, " ".join([str(e) for e in x]))
# quit()


def print_all(A, B, L, D):
    global print_count
    print "------ print {} -----".format(print_count)
    print 'A: '
    print numpy.array(A, dtype=numpy.int64)
    print 'B: '
    print numpy.array(B, dtype=numpy.int64)
    print 'L: '
    print numpy.array(L, dtype=numpy.int64)
    print 'D: '
    print numpy.array(D, dtype=numpy.int64)
    print_count += 1


if __name__ == '__main__':

    Ab = numpy.loadtxt(os.path.join(os.environ['HOME'], 'Desktop',
                                    'test_hermite.txt'))
#     print '$test_hermite = "{}";'.format(
#         " ".join([str(e) for e in Ab.ravel()]))
    x = solve(Ab[:, :-1], Ab[:, -1])
    print "The solution is x: {}".format(x)

#     offset = 0
#     if offset:
#         end = offset + 1
#     else:
#         end = 5
#     # end = 1
#     for count, (arr, x) in enumerate(zip(arrays[offset:end], xs[offset:end])):
#         print "\n\n-------- {} ----------".format(count + offset)
#         arr = arr[:-3, :-3]
#         x = x[:-2]
#         Ab = arr.T
#         G = numpy.concatenate(
#             (Ab, numpy.zeros((Ab.shape[0], 1), dtype=numpy.int64)), axis=1)
#         G[-1, -1] = 1
#         A, B, L, D = initialise_working_matrices(G)
#     #     print_all(A, B, L, D)
#         k = 3
#         i = k - 1
#         j = 3
    # Swap
    #     print "swap2($k, $m, $n): "
    #     swap_rows(k, A, B, L, D)
    #     print_all(A, B, L, D)
    # Reduce:
    #     col1, col2 = reduce_matrix(A, B, L, k, i, D)
    #     print "reduce2({k}, {i}, {m}, {n}, D): {col1}, {col2}".format(
    #         k=k, i=i, m=A.shape[0], n=A.shape[1], col1=col1, col2=col2)
    #     print_all(A, B, L, D)
    #     minus(j, A[:A.shape[1], :])
    #     print "minus(j, m, L): "
    #     print_all(A, B, L, D)
    # Hermite:
    #     hnf, unimodular_matrix, rank = lllhermite(G, m1=1, n1=1)
    #     print "lllhermite(G, {}, {}, 1, 1): {} ".format(
    #         A.shape[0], A.shape[1], rank)
    #     print "HNF:"
    #     print hnf
    #     print "Unimodular matrix:"
    #     print unimodular_matrix
    #     print arr
    # Gram:
    #     G = gram(arr)
    #     print "G:"
    #     print G
    # LCV:
    #     print "A:"
    #     print arr.T
    #     print "x:"
    #     print x
    #     print "lcv: " + str(lcasvector(arr.T, x))
    # Cholesky:
    #     PD = numpy.dot(arr, arr.T) + 1
    #     print "PD:"
    #     print PD
    #     N, D = cholesky(PD)
    #     print "cholesky(G):"
    #     print "N:"
    #     print N
    #     print "D:"
    #     print D
    # Solve:
    
    
    #     print "shortest_distance(A, m, n): " + shortest_distance(A, m, n)
    
    #     a = [-2, -1, 9, 1, 2]
    #     b = [4, 2, -5, 7, -6]
    #     c = [8, -1, 11, -1, 5]
    #     d = [-5, 1, 3, -2, 1]
    #
    #     for i in xrange(5):
    #         print "-------------- i = " + str(i) + " --------------"
    #         print "introot(" + str(a[i]) + ", " + str(b[i]) + ", " + str(c[i]) + ", " + str(d[i]) + "): {}".format(
    #             introot(abs(a[i]), abs(b[i]), c[i], d[i]))
    #         print "egcd(" + str(a[i]) + ", " + str(b[i]) + "): {}, {}, {}".format(
    #             *egcd(a[i], b[i]))
    #         print "lnearint(" + str(a[i]) + ", " + str(b[i]) + "): {}".format(
    #             lnearint(a[i], b[i]))
    #         print "ratior(" + str(a[i]) + ", " + str(b[i]) + ", " + str(c[i]) + ", " + str(d[i]) + "): {}, {}".format(
    #             *ratior(a[i], b[i], c[i], d[i]))
    #         print "multr(" + str(a[i]) + ", " + str(b[i]) + ", " + str(c[i]) + ", " + str(d[i]) + "): {}, {}".format(
    #             *multr(a[i], b[i], c[i], d[i]))
    #         print "subr(" + str(a[i]) + ", " + str(b[i]) + ", " + str(c[i]) + ", " + str(d[i]) + "): {}, {}".format(
    #             *subr(a[i], b[i], c[i], d[i]))
    #         print "addr(" + str(a[i]) + ", " + str(b[i]) + ", " + str(c[i]) + ", " + str(d[i]) + "): {}, {}".format(
    #             *addr(a[i], b[i], c[i], d[i]))
    #         print "comparer(" + str(a[i]) + ", " + str(b[i]) + ", " + str(c[i]) + ", " + str(d[i]) + "): {}".format(
    #             comparer(a[i], abs(b[i]), c[i], abs(d[i])))

