from builtins import zip
from unittest import TestCase
from diophantine import solve
from sympy import Matrix


class TestDiophantine(TestCase):

    def test_dimension_basis(self):
        """
        This test comes from the mapping of compound dimensions (b) onto a new
        set of basis dimensions (A), where each row corresponds to the 7 basic
        dimensions ('mass', 'length', 'time', 'current', 'amount',
        'temperature' and 'luminous intensity')

        The test cases correspond to:
            b_names = ['mV_per_ms', 'uF_per_cm2', 'uF_uS',
                       'K3_ms_uS4_per_cd2_cm', 'A2_mV_per_ms_uS',
                       'ms_uF3_um_per_mV', 'um3_per_C2', 'K2_per_mV_ms_ohm_uF']
            A_names = [
                ['ms', 'mV', 'mA_per_cm2', 'nA', 'mM', 'uF_per_cm2', 'um',
                 'S_per_cm2', 'uS', 'cm_ohm', 'ohm', 'degC', 'cd'],
                ['ms', 'mV', 'pA', 'mM', 'uF', 'um', 'uS', 'degC', 'cd']]
        """
        As = [
            Matrix([
                [0, 1, 0, 0, 0, -1, 0, -1, -1, 1, 1, 0, 0],
                [0, 2, -2, 0, -3, -4, 1, -4, -2, 3, 2, 0, 0],
                [1, -3, 0, 0, 0, 4, 0, 3, 3, -3, -3, 0, 0],
                [0, -1, 1, 1, 0, 2, 0, 2, 2, -2, -2, 0, 0],
                [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]]),
            Matrix([
                [0, 1, 0, 0, -1, 0, -1, 0, 0],
                [0, 2, 0, -3, -2, 1, -2, 0, 0],
                [1, -3, 0, 0, 4, 0, 3, 0, 0],
                [0, -1, 1, 0, 2, 0, 2, 0, 0],
                [0, 0, 0, 1, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 1, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 1]])]
        bs = [Matrix([1, 2, -4, -1, 0, 0, 0]),
              Matrix([-1, -4, 4, 2, 0, 0, 0]),
              Matrix([-2, -4, 7, 4, 0, 0, 0]),
              Matrix([-4, -9, 13, 8, 0, 3, -2]),
              Matrix([2, 4, -7, -1, 0, 0, 0]),
              Matrix([-4, -7, 16, 7, 0, 0, 0]),
              Matrix([0, 3, -2, -2, 0, 0, 0]),
              Matrix([-1, -2, 1, 1, 0, 2, 0])]
        num_sols = [(2, 2), (1, 1), (1, 1), (3, 1), (1, 2), (1, 1), (3, 1),
                    (1, 1)]
        for b, nsols in zip(bs, num_sols):
            for A, nsol in zip(As, nsols):
                sols = solve(A, b)
                # Check number of solutions matches reference
                self.assertEquals(len(sols), nsol,
                                  "Incorrect number of solutions found ({}), "
                                   "expected {}".format(len(sols), nsol))
                for sol in sols:
                    self.assertEqual(b, A * sol, "A * x doesn't match b")

    def test_random(self):

        As = [
            Matrix([
                [0, 3, -7, 7, -5, 4, 4, -1, -5, -9],
                [-9, -2, -1, 9, -6, 9, 1, 8, -1, 8],
                [4, 4, -4, 2, 4, 2, 5, 3, 9, 0],
                [4, 3, -5, 9, -2, 1, -7, 2, 2, 8],
                [7, 6, 5, -2, -9, -2, 0, 6, -2, -3]]),
            Matrix([
                [-5, 8, 4, -6, 7, 1, 0, 5, 8, -3],
                [1, -8, 0, 7, -4, 2, 4, -8, 5, 1],
                [-2, 0, -1, 5, -3, -2, 8, -4, -3, 8],
                [0, 5, -6, 1, 2, -3, -2, -5, -1, -9],
                [-2, 3, 1, -1, 7, -5, -9, -5, 4, -4]]),
            Matrix([
                [-5, 2, 7, -8, 3, -5, -8, -8, -1, -5],
                [-3, 5, 3, 5, 3, -8, -6, -6, -1, 5],
                [-5, -5, -5, 9, 9, 0, -2, 2, 5, -7],
                [0, 0, 1, -1, 1, -4, 9, 7, 8, 4],
                [8, -3, -1, 3, 1, 8, 6, 0, -2, -5]]),
            Matrix([
                [-5, 1, -1, 5, -3, 0, -7, 4, -9, 5],
                [-3, -6, 8, 3, 1, -7, 5, -3, 2, 3],
                [-7, 8, 3, -3, -7, 9, -5, -8, -8, 2],
                [6, -2, -3, -8, 1, -8, -4, -4, -7, 8],
                [0, -6, 1, 8, -6, -1, -1, -4, 4, 4]]),
            Matrix([
                [-5, -2, 2, -7, 3, 0, -6, -8, -1, -5],
                [-3, 0, 3, 4, 3, 4, 2, 3, -3, -9],
                [9, 0, -1, 5, 8, -2, 4, 8, -9, 1],
                [4, 3, 1, 2, 8, -6, 0, 0, 9, 3],
                [2, 2, -8, 0, 6, -2, -8, -6, 4, -9]]),
            Matrix([
                [4, -9, 4, -7, 8, -1, 1, -9, -8, -6],
                [-7, -2, 7, -4, -7, 9, 5, 4, -8, 3],
                [-2, -2, 5, 8, -5, 8, 5, -1, 3, 5],
                [-4, -4, 7, 2, -2, 2, 1, 7, -9, 2],
                [-8, -9, -4, -4, 1, 0, 2, -5, -5, 6]]),
            Matrix([
                [-5, 8, -2, -5, -1, -8, -5, -1, 5, 2],
                [-4, 3, -5, -2, -9, -8, 2, -8, 8, -1],
                [-2, 3, 0, -6, 2, 3, -1, 2, 9, -6],
                [9, -4, 1, -7, -1, 3, 2, 4, 6, 6],
                [-9, 6, -1, -8, 2, 1, 5, 5, -8, 8]]),
            Matrix([
                [9, 4, -7, -5, -4, 5, -7, 8, -5, -3],
                [7, -9, -2, -9, 8, 1, -6, -9, -3, -2],
                [-4, 4, 2, 7, -1, -5, 0, -5, -7, -9],
                [-8, 5, 6, -9, 8, 4, 7, -4, -1, 5],
                [3, -8, -6, 2, 8, -3, 9, -9, -9, -4]]),
            Matrix([
                [-9, 9, 9, 8, 2, 2, -8, 8, 4, -8],
                [-3, -7, -6, 6, -4, 7, 5, -6, 1, 1],
                [-9, 4, -2, 9, 9, 6, -5, 7, 8, 2],
                [7, -9, -5, 6, -2, 6, 6, 4, 2, 7],
                [7, 4, 9, 8, -4, -4, -9, 1, -9, 0]]),
            Matrix([
                [-2, 2, 4, 9, 3, 9, -5, -7, -3, 5],
                [-1, -2, 7, -2, 9, -2, 3, 9, -9, 0],
                [0, 0, 6, 0, -3, -9, 7, -2, 8, -4],
                [9, -2, -3, -3, 6, 6, -8, -8, 9, -6],
                [-7, -2, -5, 1, 9, 7, 3, 5, -8, -9]])]
        bs = [Matrix([-2, 0, -7, 1, 2]),
              Matrix([2, 6, -3, -6, 8]),
              Matrix([-8, 4, -8, -5, -2]),
              Matrix([-1, -9, 0, -4, 5]),
              Matrix([7, -6, -4, 1, 0]),
              Matrix([-7, 4, 7, -8, 9]),
              Matrix([-1, 0, 2, -8, 6]),
              Matrix([-9, 6, 7, -7, -7]),
              Matrix([-6, -3, -4, -8, 5]),
              Matrix([-2, -2, 6, 2, 4])]
        nsols = [1, 1, 2, 1, 1, 1, 1, 5, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 0,
                 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 2, 1,
                 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4,
                 1, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        count = 0
        for A in As:
            for b in bs:
                xs = solve(A, b)
                self.assertEqual(len(xs), nsols[count],
                                 "Incorrect number of solutions")
                for x in xs:
                    self.assertEquals(A * x, b, "A * x doesn't match b")
                count += 1


if __name__ == '__main__':
    test = TestDiophantine('test_dimension_basis')
    test.test_random()
