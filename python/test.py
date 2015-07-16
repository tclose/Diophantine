from unittest import TestCase
import os
import numpy
from nineml import units as un
from pype9.neuron import units as neuron_units
from pype9.nest import units as nest_units
from diophantine import solve
from operator import add
from sympy import Matrix


class TestDiophantine(TestCase):

    test_units = [un.mV / un.ms,
                  un.uF / un.cm ** 2,
                  un.uS * un.uF,
                  un.uS ** 4 * un.K ** 3 / (un.cd ** 2 * un.cm_per_ms),
                  un.A ** 2 * un.mV / (un.ms * un.uS),
                  un.um * un.uF ** 3 * un.ms / un.mV,
                  un.um ** 3 / un.C ** 2,
                  un.K ** 2 / (un.ohm * un.mV * un.ms * un.uF)]

    def test_solve(self):
        for unit in self.test_units:
            b = numpy.array(list(unit.dimension))
            print "Unit '{}':".format(unit.name)
            for sim_name, basis in (('NEURON', neuron_units),
                                    ('NEST', nest_units)):
                A = numpy.array([numpy.array(list(u.dimension),
                                             dtype=numpy.int64)
                                 for u in basis]).T
                solutions = solve(A, b)
                print solutions
                for solution in solutions:
                    print '  {}: '.format(sim_name) + ', '.join(
                        '{}={}'.format(u.name, d)
                        for u, d in zip(basis, solution))
                    new_dim = reduce(add, (numpy.array(list(u.dimension)) * v
                                           for u, v in zip(basis, solution)))
                    self.assertEqual(list(new_dim), list(unit.dimension),
                                     "Reconstructed dimension ({}) does not "
                                     "match original ({})".format(
                                         list(new_dim), list(unit.dimension)))


test_matrices = [
    Matrix([
        [-3, -2, -4, -3, -1, 0, -3, 0, 1, 3],
        [3, -4, 3, -1, 3, -2, -4, -2, -1, 0],
        [2, 1, 0, -2, -4, 3, 3, -4, 0, 0],
        [4, 4, 3, -4, 2, 4, 1, 0, -3, -2],
        [1, 2, 2, 1, -2, 0, 2, 0, -3, -1],
        [4, 0, -2, -1, 0, 4, 4, 2, 0, 0],
        [-4, 1, -4, 4, -4, 0, -2, 3, 4, 4]]),
    Matrix([
        [-1, 0, 0, -3, -3, -3, 4, 0, 1, 4],
        [0, -2, -2, 4, 2, -4, 0, -3, -4, 2],
        [-2, 3, 1, -4, 2, -1, 1, -4, 0, 1],
        [4, -3, -2, 2, -1, 1, -4, -2, 4, 1],
        [-3, -2, -1, -3, 0, -4, 1, -3, 3, 1],
        [-4, -1, 0, -3, 0, 0, 3, 3, -4, 0],
        [1, 1, -1, -3, 2, 2, -3, 3, 2, 2]]),
    Matrix([
        [-1, -2, 3, 0, 4, 0, -4, -3, 4, -2],
        [-3, 2, -2, 0, -4, 3, 3, 2, 0, -4],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-1, -1, 0, 0, -1, -1, 3, -3, 4, 2],
        [3, -1, 2, 0, -4, -3, -1, -1, 2, 3],
        [-4, 0, -1, 0, 4, 0, 1, -4, -2, 0],
        [4, 2, 3, 0, 0, -2, 2, -2, -4, 1]]),
    Matrix([
        [-3, 3, 4, -1, 0, -4, -1, -4, 2, -2],
        [1, 2, 3, -1, -3, 3, -3, -2, 1, -2],
        [-4, 2, 2, -2, -3, -1, -2, -4, 0, 2],
        [0, -4, -3, -3, 1, 2, 0, -3, 1, -1],
        [-1, -1, 3, 1, 1, 4, -3, -3, 0, 2],
        [0, 1, -4, 1, -3, 0, -1, 0, 1, 0],
        [0, 0, -2, -2, 4, 0, 4, 1, 2, 0]]),
    Matrix([
        [4, -3, 0, 0, 0, 3, 4, -4, 0, -3],
        [-4, 4, -3, 0, -3, -3, 2, 0, -1, -1],
        [0, 3, 4, 0, 2, -2, 2, 2, 0, 3],
        [-3, 1, 0, 0, 2, 0, 0, -3, 1, 1],
        [0, -4, -3, 0, 0, 1, -3, -1, 1, 0],
        [4, 3, 2, 0, 1, -1, 0, -2, 2, -2],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])]

xs = [
    Matrix([3, 4, 1, 0, -3, 0, 1, 4, -2]),
    Matrix([0, 4, 4, 3, 2, 4, 0, 1, 4]),
    Matrix([0, -4, 0, 1, -4, 4, 4, -2, 3]),
    Matrix([1, 2, -3, 3, -4, 1, -3, -3, -4]),
    Matrix([0, -4, -1, -2, -4, 0, 4, 3, -4])]


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
    print A
    print 'B: '
    print B
    print 'L: '
    print L
    print 'D: '
    print D
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


if __name__ == '__main__':
    test = TestDiophantine('test_solve')
    test.test_solve()
