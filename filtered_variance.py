"""Reference recursions for the deterministic filtered interval model.

Canonical distributions and sign-law moments use exact rational arithmetic.
The Gram recursion uses floating-point Euclidean norms.
"""

from collections import defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from fractions import Fraction
from math import hypot, isfinite, pi, sqrt

Number = int | Fraction
Distribution = dict[Fraction, Fraction]
C2 = (1 + 1 / sqrt(2)) / 2
GAMMA = 3 / (15 - 8 * sqrt(2))
C3_LOWER = (1 + 4 / pi) / 3
EXTREMAL_HEAVY = (1 + 3 * sqrt(2)) / 2


@dataclass(frozen=True)
class Score:
    """A score of deterministic weight on the discrete interval (issue, maturity]."""

    issue: int
    maturity: int
    weight: Number


@dataclass(frozen=True)
class LaminarOptima:
    distribution: Distribution
    variance: Fraction
    covariance_sdp: float


def _scores(items: Iterable[Score]) -> tuple[Score, ...]:
    scores = tuple(items)
    seen = set()
    for score in scores:
        if type(score.issue) is not int or type(score.maturity) is not int:
            raise TypeError("The clock uses integer endpoints.")
        if score.issue >= score.maturity:
            raise ValueError("Each interval has positive length.")
        if type(score.weight) not in (int, Fraction):
            raise TypeError("Use integer or Fraction weights for exact masses.")
        interval = score.issue, score.maturity
        if interval in seen:
            raise ValueError("Interval supports must be distinct.")
        seen.add(interval)
    return scores


def _contains(outer: Score, inner: Score) -> bool:
    return outer.issue <= inner.issue and inner.maturity <= outer.maturity


def _disjoint(left: Score, right: Score) -> bool:
    return left.maturity <= right.issue or right.maturity <= left.issue


def _length(score: Score) -> int:
    return score.maturity - score.issue


def _convolve(left: Distribution, right: Distribution) -> Distribution:
    combined: Distribution = defaultdict(Fraction)
    for a, probability_a in left.items():
        for b, probability_b in right.items():
            combined[a + b] += probability_a * probability_b
    return dict(combined)


def _extend_at_root(children: Distribution, magnitude: Fraction) -> Distribution:
    law: Distribution = defaultdict(Fraction)
    for value, mass in children.items():
        if value == 0:
            law[-magnitude] += mass / 2
            law[magnitude] += mass / 2
        else:
            law[value + (magnitude if value > 0 else -magnitude)] += mass
    return dict(law)


def laminar_optima(items: Iterable[Score]) -> LaminarOptima:
    """Return the attained forest law, its exact variance, and a numerical SDP value.

    The support size may grow exponentially. An empty forest has both optima zero.
    For signed input, deterministic sign flips reduce the weights to magnitudes.
    """
    scores = _scores(items)
    children: dict[Score, list[Score]] = {score: [] for score in scores}
    maximal = []
    for i, score in enumerate(scores):
        parents = []
        for j, candidate in enumerate(scores):
            if i == j or _disjoint(score, candidate):
                continue
            if _contains(candidate, score):
                parents.append(candidate)
            elif not _contains(score, candidate):
                raise ValueError("Crossing intervals are outside the laminar recursion.")
        if parents:
            children[min(parents, key=_length)].append(score)
        else:
            maximal.append(score)

    subtree_law: dict[Score, Distribution] = {}
    subtree_norm: dict[Score, float] = {}
    zero = {Fraction(0): Fraction(1)}
    for score in sorted(scores, key=_length):
        law = zero
        for child in children[score]:
            law = _convolve(law, subtree_law[child])
        magnitude = abs(Fraction(score.weight))
        subtree_law[score] = _extend_at_root(law, magnitude)
        subtree_norm[score] = float(magnitude) + hypot(
            *(subtree_norm[child] for child in children[score])
        )

    law = zero
    for component in maximal:
        law = _convolve(law, subtree_law[component])
    variance = sum((value**2 * mass for value, mass in law.items()), Fraction(0))
    sdp = sum((subtree_norm[component] ** 2 for component in maximal), 0.0)
    return LaminarOptima(law, variance, sdp)


def sign_law_variance(
    items: Iterable[Score], masses: Mapping[tuple[int, ...], Number]
) -> Fraction:
    """Verify all matured-history LP equations and return the exact objective.

    Unlike the laminar recursion, this checker accepts crossing intervals.
    It validates a supplied law over sign atoms.
    """
    scores = _scores(items)
    histories = tuple(
        tuple(j for j, past in enumerate(scores) if past.maturity <= score.issue)
        for score in scores
    )
    mass_total = Fraction(0)
    objective = Fraction(0)
    balances: dict[tuple[int, tuple[int, ...]], Fraction] = defaultdict(Fraction)
    for signs, raw_mass in masses.items():
        if len(signs) != len(scores) or any(type(z) is not int or z not in (-1, 1) for z in signs):
            raise ValueError("Each atom must contain one sign per score.")
        if type(raw_mass) not in (int, Fraction):
            raise TypeError("Use integer or Fraction masses.")
        mass = Fraction(raw_mass)
        if mass < 0:
            raise ValueError("Atom probabilities must be nonnegative.")
        mass_total += mass
        total = sum((score.weight * z for score, z in zip(scores, signs)), Fraction(0))
        objective += mass * total**2
        for i, history in enumerate(histories):
            key = i, tuple(signs[j] for j in history)
            balances[key] += mass * signs[i]
    if mass_total != 1:
        raise ValueError("The atom probabilities must sum to one.")
    if any(balances.values()):
        raise ValueError("A matured-history conditional centering equation fails.")
    return objective


def paired_slice(heavy: float, *, like_with_like: bool = False) -> tuple[float, float]:
    """Return (optimal outer weight, minimum ratio) for (heavy, heavy, 1, 1).

    This is a symmetric slice of the seven-weight problem. The two matching
    options have the same stochastic moments. Each is optimized over r >= 0.
    """
    if not isfinite(heavy) or heavy < 1:
        raise ValueError("The heavy coefficient must be finite and at least one.")
    absolute_moment = heavy + 0.5
    second_moment = 2 * (heavy**2 + 1)
    if like_with_like:
        gram_length = sqrt(second_moment)
    else:
        gram_length = sqrt(2) * (heavy + sqrt(2) - 1)
    outer_weight = max(0.0, (second_moment - absolute_moment * gram_length)
                       / (gram_length - absolute_moment))
    ratio = (outer_weight**2 + 2 * outer_weight * absolute_moment + second_moment)
    ratio /= (outer_weight + gram_length)**2
    return outer_weight, ratio
