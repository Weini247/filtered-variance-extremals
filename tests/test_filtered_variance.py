from collections import Counter
from fractions import Fraction
from itertools import product
import unittest

from filtered_variance import C2, EXTREMAL_HEAVY, GAMMA, Score, laminar_optima, paired_slice, sign_law_variance


class FilteredVarianceTests(unittest.TestCase):
    def test_disjoint_scores_and_shared_issue_maturity(self):
        answer = laminar_optima([Score(0, 1, -2), Score(1, 2, 3)])
        self.assertEqual(answer.variance, 13)
        self.assertEqual(answer.covariance_sdp, 13)
        self.assertEqual(answer.distribution, {z: Fraction(1, 4) for z in (-5, -1, 1, 5)})

    def test_root_uses_fair_coin_at_zero(self):
        scores = [Score(0, 2, 1), Score(0, 1, 1), Score(1, 2, 1)]
        answer = laminar_optima(scores)
        self.assertEqual(answer.distribution, {z: Fraction(1, 4) for z in (-3, -1, 1, 3)})
        self.assertEqual(answer.variance, 5)
        self.assertAlmostEqual(answer.covariance_sdp, 3 + 2**1.5)

    def test_seven_unit_weights(self):
        scores = [Score(0, 6, 1), Score(0, 3, 1), Score(0, 1, 1), Score(1, 3, 1),
                  Score(3, 6, 1), Score(3, 4, 1), Score(4, 6, 1)]
        answer = laminar_optima(scores)
        self.assertEqual(answer.variance, 16)
        self.assertAlmostEqual(answer.covariance_sdp, 11 + 6 * 2**0.5)
        self.assertLess(float(answer.variance) / answer.covariance_sdp, C2)

    def test_depth_four_fixture(self):
        items = (Score(0, 5, 1), Score(0, 3, 2), Score(0, 2, 1), Score(0, 1, 1),
                 Score(1, 2, 1), Score(2, 3, 2), Score(3, 4, 0))
        answer = laminar_optima(items)
        self.assertEqual(answer.variance, 33)
        self.assertEqual(answer.distribution, {-8: Fraction(1, 8), -6: Fraction(1, 8),
                         -4: Fraction(1, 4), 4: Fraction(1, 4), 6: Fraction(1, 8),
                         8: Fraction(1, 8)})

    def test_fractional_and_zero_weights(self):
        answer = laminar_optima([Score(0, 2, Fraction(1, 2)), Score(0, 1, 0)])
        self.assertEqual(answer.variance, Fraction(1, 4))
        self.assertEqual(laminar_optima([]).distribution, {0: 1})
        self.assertEqual(laminar_optima([Score(0, 1, 0)]).distribution, {0: 1})

    def test_bad_interval_geometry(self):
        cases = ([Score(0, 3, 1), Score(1, 4, 1)],
                 [Score(2, 2, 1)],
                 [Score(0, 2, 1), Score(0, 2, 0)])
        for scores in cases:
            with self.subTest(scores=scores), self.assertRaises(ValueError):
                laminar_optima(scores)

    def test_feasible_sign_law_matches_recursion(self):
        scores = [Score(0, 2, -1), Score(0, 1, 1), Score(1, 2, 1)]
        counts = Counter()
        for first, second, tie in product((-1, 1), repeat=3):
            magnitude_sign = first if first == second else tie
            counts[(-magnitude_sign, first, second)] += 1
        law = {atom: Fraction(count, 8) for atom, count in counts.items()}
        self.assertEqual(sign_law_variance(scores, law), laminar_optima(scores).variance)

    def test_pairwise_orthogonal_parity_is_not_feasible(self):
        scores = [Score(i, i + 1, 1) for i in range(3)]
        law = {(r, s, r*s): Fraction(1, 4) for r, s in product((-1, 1), repeat=2)}
        for i, j in ((0, 1), (0, 2), (1, 2)):
            self.assertEqual(sum(atom[i]*atom[j]*mass for atom, mass in law.items()), 0)
        with self.assertRaisesRegex(ValueError, "centering"):
            sign_law_variance(scores, law)

    def test_crossing_sign_law_is_part_of_the_lp_domain(self):
        scores = [Score(0, 2, 1), Score(1, 3, 1)]
        self.assertEqual(sign_law_variance(scores, {(1, 1): Fraction(1, 2),
                                                    (-1, -1): Fraction(1, 2)}), 4)

    def test_figure_slice_at_the_analytic_witness(self):
        outer, ratio = paired_slice(EXTREMAL_HEAVY)
        self.assertAlmostEqual(outer, 2)
        self.assertAlmostEqual(ratio, GAMMA)
        self.assertGreater(paired_slice(EXTREMAL_HEAVY, like_with_like=True)[1], ratio)
        for mixed, same in zip(paired_slice(1), paired_slice(1, like_with_like=True)):
            self.assertAlmostEqual(mixed, same)

    def test_fractional_lp_mass_checks(self):
        score = [Score(0, 1, 1)]
        for masses in ({(1,): Fraction(1, 2)}, {(1,): Fraction(-1), (-1,): 2},
                       {(0,): 1}):
            with self.subTest(masses=masses), self.assertRaises(ValueError):
                sign_law_variance(score, masses)


if __name__ == "__main__":
    unittest.main()
