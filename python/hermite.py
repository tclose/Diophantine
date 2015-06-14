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
"""
hermite.py

Input: integer mxn matrix A, nonzero, at least two rows+
Output: small unimodular matrix B and HNF(A), such that BA=HNF(A)+
The Havas, Majewski, Matthews LLL method is used+
We usually take alpha=m1/n1, with (m1,n1)=(1,1) to get best results+
"""

# global col1
# global col2
# global n + 1
# global B
# global L
# global A
# global D
# global hnf
# global unimodular_matrix
# global rank
from copy import deepcopy
from fractions import gcd
import numpy


def axb(Ab, m, n, m1, n1):
    """
    A is m x n, b is m x 1, solving AX=b, X is n x 1+
    Ab is the (n+1) x m transposed augmented matrix. G=[A^t|0] [b^t]1]
    """
    # global hnf
    # global unimodular_matrix
    # global rank
#     G = deepcopy(Ab)
    for i in xrange(m + 1):
        for j in xrange(n):
            G[i][j] = Ab[i][j]
    for i in xrange(m):
        G[i][n + 1] = 0
    G[m + 1][n + 1] = 1
#    print "G="
#    printmat1(G,m + 1,n + 1)
#    print "<br>\n"
    hnf, unimodular_matrix, rank = lllhermite(G, m + 1, n + 1, m1, n1)
#    print "HNF(G)="
#    printmat1(hnf,m + 1,n + 1)
#    print "<br>\n"
#    print "P ="
#    printmat1(unimodular_matrix,m + 1,m + 1)
#    print "is a unimodular matrix such that PG = HNF(G)"
#    print "<br>\n"
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
        # print "<img align=\"middle\" src=\"../jpgs/matrixP.png\"><br>\n"
        y = -unimodular_matrix[rank, :]
        # print "AX=B has a solution: Y = "
        # print[y,m]
        # print "<br>\n"
        nullity = m + 1 - rank
        if nullity == 0:
            print "AX=B has a unique solution in integers<br>\n"
            return
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
    shortest_distance_axb(basis, lim + 1, m)
    return


def lllhermite(G, m, n, m1, n1):
    """  G is a nonzero matrix with at least two rows.  """
    m = G.shape[0]
    n = G.shape[1]
    B = numpy.eye(m)
    L = numpy.zeros((m, m))  # Lower triangular matrix
    D = numpy.ones(m + 1)
    for i in xrange(m):
        for j in xrange(n):
            A[i][j] = G[i][j]
    flag = flagcol(A, m, n)
    if flag == 1:
        B[m][m] = -1
        for j in xrange(n):
            A[m][j] = -A[m][j]
    k = 2
    while k <= m:
        reduce2(k, k - 1, m, n, D)
        minim = min(col2, n)
        temp1 = D[k - 2] * D[k]
        temp2 = L[k][k - 1] * L[k][k - 1]
        temp3 = temp1 + temp2
        u = n1 * temp3
        temp1 = D[k - 1] * D[k - 1]
        v = m1 * temp1
        if col1 <= minim or (col1 == col2 and col1 == n + 1 and u < v):
            swap2(k, m, n)
            if k > 2:
                k = k - 1
        else:
            for i in xrange(k - 2 - 1, 0, -1):
                reduce2(k, i, m, n, D)
            k = k + 1
    for i in xrange(m):
        for j in xrange(n):
            hnf[i][j] = A[i][j]
    for i in xrange(m):
        for j in xrange(m):
            unimodular_matrix[i][j] = B[i][j]
    for i in xrange(m - 1, 0, -1):
        test = zero_row_test(A, n, i)
        if test == 0:
            break
    rank = m - i
    for i in xrange(m):
        for j in xrange(n):
            k = m + 1 - i
            hnf[i][j] = A[k][j]
    for i in xrange(m):
        for j in xrange(m):
            k = m + 1 - i
            unimodular_matrix[i][j] = B[k][j]
    return hnf, unimodular_matrix, rank


def flagcol(A, m, n):
    """
    returns 0 if the first nonzero column j of A contains more than one nonzero
    entry, or contains only one nonzero entry and which is positive+ returns 1
    if the first nonzero column j of A contains only one nonzero entry, which
    is negative+ This assumes A is a nonzero matrix with at least two rows+
    """
    flag = 0
    for j in xrange(n):
        for i in xrange(m):
            # found the first column with a nonzero elt, which is in row i
            if A[i][j] != 0:
                flag = 1
                break
        if flag == 1:
            break
    for k in xrange(i, m):
        if A[k][j] != 0:
            return 0
    if A[i][j] > 0:  # A[i][j] is the only elt in column j and is positive
        return 0
    else:  # A[i][j] is the only elt in column j and is negative
        return 1


def reduce2(k, i, m, n, D):
    # global col1
    # global col2
    # global n + 1
    # global B
    # global L
    # global A
    col1 = n + 1
    for j in xrange(n):
        if A[i][j] != 0:
            col1 = j
        if A[i][col1] < 0:
            minus(i, m, L)
            for jj in xrange(n):
                A[i][jj] = -A[i][jj]
            for jj in xrange(m):
                B[i][jj] = -B[i][jj]
        break
    col2 = n + 1
    for j in xrange(n):
        if A[k][j] != 0:
            col2 = j
            break
    if col1 <= n:
        q = A[k][col1] // A[i][col1]
    else:
        t = abs(L[k][i])
        t = 2 * t
        if t > D[i]:
            q = lnearint(L[k][i], D[i])
        else:
            q = 0
    if q != 0:
        for j in xrange(n):
            temp = q * A[i][j]
            A[k][j] = A[k][j] - temp
        for j in xrange(m):
            temp = q * B[i][j]
            B[k][j] = B[k][j] - temp
        temp = q * D[i]
        L[k][i] = L[k][i] - temp
        for j in xrange(i - 1):
            temp = q * L[i][j]
            L[k][j] = L[k][j] - temp


def minus(j, m, L):
    for r in xrange(1, m):
        for s in xrange(r - 1):
            if r == j or s == j:
                L[r][s] = -L[r][s]


def swap2(k, m, n):
    # global B
    # global L
    # global A
    # global D
    # print "Row k <. Row k - 1<br>\n"
    for j in xrange(n):
        temp = A[k][j]
        A[k][j] = A[k - 1][j]
        A[k - 1][j] = temp
    for j in xrange(m):
        temp = B[k][j]
        B[k][j] = B[k - 1][j]
        B[k - 1][j] = temp
    for j in xrange(k - 2):
        temp = L[k][j]
        L[k][j] = L[k - 1][j]
        L[k - 1][j] = temp
    for i in xrange(k, m):
        temp1 = L[i][k - 1] * D[k]
        temp2 = L[i][k] * L[k][k - 1]
        t = temp1 - temp2
        temp1 = L[i][k - 1] * L[k][k - 1]
        temp2 = L[i][k] * D[k - 2]
        temp3 = temp1 + temp2
        L[i][k - 1] = temp3 / D[k - 1]
        L[i][k] = t / D[k - 1]
    temp1 = D[k - 2] * D[k]
    temp2 = L[k][k - 1] * L[k][k - 1]
    t = temp1 + temp2
    D[k - 1] = t / D[k - 1]
    return


def zero_row_test(matrix, n, i):
    """
    This tests the i-th row of matrix to see if there is a nonzero
    entry. If there is one and the first occurs in column j, then j
    is returned. Otherwise 0 is returned
    """
    for j in xrange(n):
        if matrix[i][j] != 0:
            return j
    return 0


def shortest_distance_axb(AXB):
    # global choleskynum
    # global choleskyden
    # global multnum
    # global multden
    # global addnum
    # global addden
    # global subnum
    # global subden
    # global rationum
    # global ratioden
    # global lcv
    count = 0
    m = AXB.shape[0]
    n = AXB.shape[1]
    m -= 1  # Not sure about this
    AA = AXB[:-1, :]  # AA consists of the first m-1 rows of A
    G = gram(AXB)
    Qn, Qd = cholesky(G)
    m -= 1
    for i in xrange(m):  # the N vector
        Nn[i] = Qn[i][m + 1]
        Nd[i] = Qd[i][m + 1]
    Cn = 0
    Cd = 1
    for i in xrange(m):
        n, d = multr(Nn[i], Nd[i], Nn[i], Nd[i])
        n, d = multr(n, d, Qn[i][i], Qd[i][i])
        Cn, Cd = addr(Cn, Cd, n, d)
    i = m
    Tn[m] = Cn
    Td[m] = Cd
    Un[m] = 0
    Ud[m] = 1
    while 1:
        Zn, Zd = ratior(Tn[i], Td[i], Qn[i][i], Qd[i][i])
        n, d = subr(Nn[i], Nd[i], Un[i], Ud[i])
        UB[i] = introot(Zn, Zd, n, d)
        n, d = subr(Un[i], Ud[i], Nn[i], Nd[i])
        x[i] = -introot(Zn, Zd, n, d) - 1
        while True:
            x[i] += 1
            if le(x[i], UB[i]):
                if i == 1:
                    lcv = lcasvector(AA, x)
                    count = count + 1
                    lcva[count] = lcv
                    coord[count] = x
                    for k in xrange(n):
                        temp = A[m + 1][k]
                        mulitpliers[count][k] = temp - lcv[k]
                    l = mulitpliers[count] ** 2  # lengthsquared(mulitpliers[count], n)
                    mulitpliers[count][n + 1] = l
                    continue
                else:
                    i = i - 1
                    # now update U[i]
                    Un[i], Ud[i] = 0, 1
                    for j in xrange(i, m):
                        n, d = multr(Qn[i][j], Qd[i][j], x[j], 1)
                        Un[i], Ud[i] = addr(Un[i], Ud[i], n, d)
                    # now update T[i]
                    n, d = addr(x[i + 1], 1, Un[i + 1], Ud[i + 1])
                    n, d = subr(n, d, Nn[i + 1], Nd[i + 1])
                    n, d = multr(n, d, n, d)
                    n, d = multr(Qn[i + 1][i + 1], Qd[i + 1][i + 1], n, d)
                    Tn[i], Td[i] = subr(Tn[i + 1], Td[i + 1], n, d)
                    break
            else:
                i = i + 1
                if i > m:
                    return mulitpliers
                continue


def cholesky(A):
    """
    # A is positive definite mxm
    """
    assert A.dim == 2 and A.shape[0] == A.shape[1]
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
    assert m == A.shape[1]
    B = numpy.empty(m)
    for i in xrange(m):
        for j in xrange(m):
            B[i][j] = A[i].dot(A[j])  # dotproduct(A[i], A[j], n)
    return B


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
    x = numpy.sqrt(x)
    answer = x + y
    n, d = subr(c, d, y, 1)
    n, d = subr(1, 1, n, d)
    n, d = addr(x, 1, n, d)
    n, d = multr(n, d, n, d)
    t = comparer(n, d, a, b)
    if t <= 0:
        answer = answer + 1
    return answer


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
    g = gcd(r, s)
    if s < 0:
        g = -g
    return r / g, s / g


def multr(a, b, c, d):
    # returns (a/b)(c/d)
    r = a * c
    s = b * d
    g = gcd(r, s)
    return r / g, s / g


def subr(a, b, c, d):
    t = a * d - b * c
    u = b * d
    g = gcd(t, u)
    return t / g, u / g


def addr(a, b, c, d):
    t = a * d + b * c
    u = b * d
    g = gcd(t, u)
    return t / g, u / g


def comparer(a, b, c, d):
    """Assumes b>0 and d>0.  Returns -1, 0 or 1 according as a/b <,=,> c/d+ """
    assert b > 0 and d > 0
    return numpy.sign(a * d - b * c)


def lcasvector(A, X):
    """lcv[j]=X[1]A[1][j]=...+X[m]A[m][j], 1 <= j <= n+"""
    # global lcv
    n = A.shape[0]
    lcv = numpy.empty(n)
    for j in xrange(n):
        lcv[j] = X.dot(A[:, j])
    return lcv
