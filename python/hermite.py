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
from copy import deepcopy
from fractions import gcd
import numpy
from nineml import units as un
global print_count
print_count = 0


def diophantine(compound, reference_dims, m1=1, n1=1):
    """
    Finds the minimal combination of reference dimensions to make the compound
    dimension
    """
    A = numpy.array([list(d) for d in reference_dims])
    b = numpy.array([list(compound)])
    m, n = A.shape
    Ab = numpy.concatenate((A, b))
    G = numpy.concatenate((Ab, numpy.zeros((m, 1))), axis=1)
    G[-1, -1] = 1
    # A is m x n, b is m x 1, solving AX=b, X is n x 1+
    # Ab is the (n+1) x m transposed augmented matrix. G=[A^t|0] [b^t]1]
    hnf, unimodular_matrix, rank = lllhermite(G, m1=m1, n1=n1)
    flag = 0
    for i in xrange(rank - 1):
        if hnf[i][n + 1] != 0:
            flag = 1
            break
    flag1 = 0
    for j in xrange(n):
        if hnf[rank][j] != 0:
            flag1 = 1
            break
    if flag == 0 and hnf[rank][n + 1] == 1 and flag1 == 0:
        y = -unimodular_matrix[rank, :]
        nullity = m + 1 - rank
        if nullity == 0:
            print "AX=B has a unique solution in integers<br>\n"
            return y
        else:
            lim = m + 1 - rank
            basis = []
            for i in xrange(lim):
                basis.append(unimodular_matrix[rank + i, :])
            basis = numpy.array(basis)
    else:
        raise Exception("AX=B has no solution in integers<br>\n")
    # joining basis and y
    basis = numpy.concatenate(basis, y)
    return shortest_distance_axb(basis, lim + 1, m)


def lllhermite(G, m1=1, n1=1):
    """
    Input: integer mxn matrix A, nonzero, at least two rows+
    Output: small unimodular matrix B and HNF(A), such that BA=HNF(A)+
    The Havas, Majewski, Matthews LLL method is used+
    We usually take alpha=m1/n1, with (m1,n1)=(1,1) to get best results+
    """
    global print_count
    m = G.shape[0]
    n = G.shape[1]
    A, B, L, D = initialise_working_matrices(G)
    if first_nonzero_is_negative(A):
        B[m, m] = -1
        A[m, :] *= -1
    k = 1
    while k < m:
        print "k={k}, m={m}".format(k=k, m=m)
        col1, col2 = reduce_matrix(A, B, L, k, k - 1, D)
        print "col1={col1}, col2={col2}".format(col1=col1, col2=col2)
        print_all(A, B, L, D)
        u = n1 * (D[k - 1] * D[k + 1] + L[k, k - 1] * L[k, k - 1])
        v = m1 * D[k] * D[k]
        print "u={u}, v={v}".format(u=u, v=v)
        if col1 <= min(col2, n) or (col1 == col2 and col1 == n + 1 and u < v):
            swap_rows(k, A, B, L, D)
            print_all(A, B, L, D)
            if k > 1:
                k = k - 1
                print "col1 <= minim && k > 1"
            else:
                print "col1 <= minim"
        else:
            for i in reversed(xrange(k - 1)):
                reduce_matrix(A, B, L, k, i, D)
            k = k + 1
            print "col1 > minim"
    hnf = deepcopy(A)
    unimodular_matrix = deepcopy(B)
    for i in xrange(m - 1, 0, -1):
        test = zero_row_test(A, i)
        if test == 0:
            break
    rank = m - i
    for i in xrange(m):
        for j in xrange(n):
            hnf[i, j] = A[m - i, j]
    for i in xrange(m):
        for j in xrange(m):
            unimodular_matrix[i, j] = B[m - i, j]
    return hnf, unimodular_matrix, rank


def initialise_working_matrices(G):
    """  G is a nonzero matrix with at least two rows.  """
    B = numpy.eye(G.shape[0], dtype=int)
    L = numpy.zeros((G.shape[0], G.shape[0]), dtype=int)  # Lower triang matrix
    D = numpy.ones(G.shape[0] + 1, dtype=int)
    A = numpy.array(G, dtype=int)
    return A, B, L, D


def first_nonzero_is_negative(A):
    """
    returns 0 if the first nonzero column j of A contains more than one nonzero
    entry, or contains only one nonzero entry and which is positive+ returns 1
    if the first nonzero column j of A contains only one nonzero entry, which
    is negative+ This assumes A is a nonzero matrix with at least two rows+
    """
    nonzero_columns = numpy.nonzero(numpy.sum(A, axis=0) != 0)[0]
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
        col1 = A.shape[1] - 1
    nonzero_k_elems = numpy.nonzero(A[k])[0]
    if len(nonzero_k_elems):
        col2 = nonzero_k_elems[0]
    else:
        col2 = A.shape[1] - 1
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
    print_all(A, B, L, D)
    t = L[(k + 1):, k - 1] * D[k + 1] - L[(k + 1):, k] * L[k, k - 1]
    L[(k + 1):, k - 1] = (L[(k + 1):, k - 1] * L[k, k - 1] +
                          L[(k + 1):, k] * D[k - 1]) / D[k]
    L[(k + 1):, k] = t / D[k]
    print_all(A, B, L, D)
    t = D[k - 1] * D[k + 1] + L[k, k - 1] * L[k, k - 1]
    D[k] = t / D[k]
    return


def zero_row_test(matrix, i):
    """
    This tests the i-th row of matrix to see if there is a nonzero
    entry. If there is one and the first occurs in column j, then j
    is returned. Otherwise 0 is returned
    """
    nonzero_elems = numpy.nonzero(matrix[i])[0]
    if len(nonzero_elems):
        return nonzero_elems[0]
    else:
        return -1


def shortest_distance_axb(A):
    m = A.shape[0]
    n = A.shape[1]
    m -= 1  # Not sure about this
    G = gram(A)
    Qn, Qd = cholesky(G)
    m -= 1
    Nn = Qn[:, m + 1]
    Nd = Qd[:, m + 1]
    Cn = 0
    Cd = 1
    for i in xrange(m):
        n, d = multr(Nn[i], Nd[i], Nn[i], Nd[i])
        n, d = multr(n, d, Qn[i][i], Qd[i][i])
        Cn, Cd = addr(Cn, Cd, n, d)
    i = m
    Tn = Cn
    Td = Cd
    Un = 0
    Ud = 1
    multipliers = []  # List to hold multipliers
    xs = []  # List to hold pas values of x
    while 1:
        # Calculate UB
        Zn, Zd = ratior(Tn, Td, Qn, Qd)
        n, d = subr(Nn, Nd, Un, Ud)
        UB = introot(Zn, Zd, n, d)
        # Calculate x
        n, d = subr(Un, Ud, Nn, Nd)
        x = -introot(Zn, Zd, n, d) - 1
        while True:
            x += 1
            if x <= UB:
                if i == 1:
                    lcv = lcasvector(A[:-1, :], x)
                    multiplier = A[m + 1, :n] - lcv
#                   lengthsquared(mulitpliers[count], n)
                    l = multiplier ** 2
                    multiplier[n + 1] = l
                    multipliers.append(multiplier)
                    continue
                else:
                    # Save x in list of previous xs
                    xs.append(x)
                    # now update U
                    prev_Un = Un
                    prev_Ud = Ud
                    Un, Ud = 0, 1
                    for j in xrange(i, m):
                        # Loops from back of xs
                        n, d = multr(Qn[j], Qd[j], xs[i - j], 1)
                        Un, Ud = addr(Un, Ud, n, d)
                    # now update T
                    n, d = addr(x, 1, prev_Un, prev_Ud)
                    n, d = subr(n, d, Nn[i], Nd[i])
                    n, d = multr(n, d, n, d)
                    n, d = multr(Qn[i][i], Qd[i][i], n, d)
                    Tn, Td = subr(Tn, Td, n, d)
                    i = i - 1
                    break
            else:
                i = i + 1
                if i > m:
                    return multipliers
                continue


def cholesky(A):
    """
    # A is positive definite mxm
    """
    assert A.ndim == 2 and A.shape[0] == A.shape[1]
    assert numpy.all(numpy.linalg.eigvals(A) > 0)
    m = A.shape[0]
    N = deepcopy(A)
    D = numpy.ones(A.shape)
    for i in xrange(1, m):
        for j in xrange(i, m):
            N[j][i] = N[i][j]
            D[j][i] = D[i][j]
            N[i][j], D[i][j] = ratior(N[i][j], D[i][j], N[i][i], D[i][i])
    for k in xrange(i, m):
        for l in xrange(k - 1, m):
            n, d = multr(N[k][i], D[k][i], N[i][l], D[i][l])
            N[k][l], D[k][l] = subr(N[k][l], D[k][l], n, d)
    return N, D


def gram(A):
    """
    Need to check for row and column operations
    """
    m = A.shape[0]
    B = numpy.empty((m, m))
    for i in xrange(m):
        for j in xrange(m):
            B[i][j] = A[i].dot(A[j])  # dotproduct(A[i], A[j], n)
    return numpy.array(B, dtype=int)


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
    lcv = numpy.empty(n)
    for j in xrange(n):
        lcv[j] = x.dot(A[:, j])
    return numpy.array(lcv, dtype=int)


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
        [-4, 1, -4, 4, -4, 0, -2, 3, 4, 4]], dtype=int),
    numpy.array([
        [-1, 0, 0, -3, -3, -3, 4, 0, 1, 4],
        [0, -2, -2, 4, 2, -4, 0, -3, -4, 2],
        [-2, 3, 1, -4, 2, -1, 1, -4, 0, 1],
        [4, -3, -2, 2, -1, 1, -4, -2, 4, 1],
        [-3, -2, -1, -3, 0, -4, 1, -3, 3, 1],
        [-4, -1, 0, -3, 0, 0, 3, 3, -4, 0],
        [1, 1, -1, -3, 2, 2, -3, 3, 2, 2]], dtype=int),
    numpy.array([
        [-1, -2, 3, 0, 4, 0, -4, -3, 4, -2],
        [-3, 2, -2, 0, -4, 3, 3, 2, 0, -4],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, -1, 0, 0, -1, -1, 3, -3, 4, 2],
        [3, -1, 2, 0, -4, -3, -1, -1, 2, 3],
        [-4, 0, -1, 0, 4, 0, 1, -4, -2, 0],
        [4, 2, 3, 0, 0, -2, 2, -2, -4, 1]], dtype=int),
    numpy.array([
        [-3, 3, 4, -1, 0, -4, -1, -4, 2, -2],
        [1, 2, 3, -1, -3, 3, -3, -2, 1, -2],
        [-4, 2, 2, -2, -3, -1, -2, -4, 0, 2],
        [0, -4, -3, -3, 1, 2, 0, -3, 1, -1],
        [-1, -1, 3, 1, 1, 4, -3, -3, 0, 2],
        [0, 1, -4, 1, -3, 0, -1, 0, 1, 0],
        [0, 0, -2, -2, 4, 0, 4, 1, 2, 0]], dtype=int),
    numpy.array([
        [4, -3, 0, 0, 0, 3, 4, -4, 0, -3],
        [-4, 4, -3, 0, -3, -3, 2, 0, -1, -1],
        [0, 3, 4, 0, 2, -2, 2, 2, 0, 3],
        [-3, 1, 0, 0, 2, 0, 0, -3, 1, 1],
        [0, -4, -3, 0, 0, 1, -3, -1, 1, 0],
        [4, 3, 2, 0, 1, -1, 0, -2, 2, -2],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]], dtype=int)]

xs = [
    numpy.array([3, 4, 1, 0, -3, 0, 1, 4, -2], dtype=int),
    numpy.array([0, 4, 4, 3, 2, 4, 0, 1, 4], dtype=int),
    numpy.array([0, -4, 0, 1, -4, 4, 4, -2, 3], dtype=int),
    numpy.array([1, 2, -3, 3, -4, 1, -3, -3, -4], dtype=int),
    numpy.array([0, -4, -1, -2, -4, 0, 4, 3, -4], dtype=int)]


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
    print numpy.array(A, dtype=int)
    print 'B: '
    print numpy.array(B, dtype=int)
    print 'L: '
    print numpy.array(L, dtype=int)
    print 'D: '
    print numpy.array(D, dtype=int)
    print_count += 1


offset = 0
if offset:
    end = offset + 1
else:
    end = 5
# end = 1
for count, (arr, x) in enumerate(zip(arrays[offset:end], xs[offset:end])):
    print "\n\n-------- {} ----------".format(count + offset)
    Ab = arr.T
    G = numpy.concatenate((Ab, numpy.zeros((Ab.shape[0], 1))), axis=1)
    G[-1, -1] = 1
    A, B, L, D = initialise_working_matrices(G)
#     print_all(A, B, L, D)
    k = 3
    i = k - 1
    j = 3
#     print "swap2($k, $m, $n): "
#     swap_rows(k, A, B, L, D)
#     print_all(A, B, L, D)
#     col1, col2 = reduce_matrix(A, B, L, k, i, D)
#     print "reduce2({k}, {i}, {m}, {n}, D): {col1}, {col2}".format(
#         k=k, i=i, m=A.shape[0], n=A.shape[1], col1=col1, col2=col2)
#     print_all(A, B, L, D)
#     minus(j, A[:A.shape[1], :])
#     print "minus(j, m, L): "
#     print_all(A, B, L, D)
#     print "zero_row_test(matrix, n, i): {}".format(zero_row_test(A, k))
#     X = gram(A)
#     print "gram(A, m, n): "
#     print X
#     print "lcasvector(A, X, m, nplus1): {}".format(lcasvector(A[:-1, :-1], x))
    hnf, unimodular_matrix, rank = lllhermite(G, m1=1, n1=1)
    print "lllhermite(G, $mplus1, $nplus1, $m1, $n1): " + str(rank)
    print "HNF:"
    print hnf
    print "Unimodular matrix:"
    print unimodular_matrix
#     N, D = cholesky(X)
#     print "cholesky(X, mplus1): "
#     print "Cholesky Num:"
#     print N
#     print "Cholesky Den:"
#     print D
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

