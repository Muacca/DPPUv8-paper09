"""Integer character projections on spin spherical space forms."""
from collections import Counter
from fractions import Fraction
from functools import cache
from math import comb, gcd


@cache
def divisors(n):
    return tuple(x for x in range(1, n + 1) if n % x == 0)


@cache
def mobius(n):
    return 1 if n == 1 else -sum(mobius(x) for x in divisors(n) if x != n)


@cache
def ramanujan(order, weight):
    return sum(d * mobius(order // d) for d in divisors(gcd(order, weight)))


@cache
def class_character(n, order):
    phi = sum(gcd(i, order) == 1 for i in range(1, order + 1))
    return Fraction(sum(ramanujan(order, u) for u in range(-n, n + 1, 2)), phi)


def binary_dimension(order, classes, left, right):
    total = sum(count * class_character(left, m) for m, count in classes.items())
    value = total * (right + 1) / order
    if value.denominator != 1 or value < 0:
        raise ArithmeticError(f"Invalid invariant dimension: {value}")
    return value.numerator


def lens_dimension(p, a, b, left, right):
    residues = Counter(a * u % p for u in range(-left, left + 1, 2))
    return sum(residues[-b * v % p] for v in range(-right, right + 1, 2))


def lifts():
    """Each callback counts invariants in V_left tensor V_right."""
    from functools import partial
    cyclic = [
        ("S3", 1, 0, 0), ("RP3_even", 2, 1, 0), ("RP3_odd", 2, 0, 1),
        ("C3_left", 3, 1, 0), ("C3_right", 3, 0, 1),
        ("L5", 5, 2, 1), ("L7_23", 7, 3, 2),
    ]
    for name, p, a, b in cyclic:
        yield name, p, partial(lens_dimension, p, a, b)
    tet = {1: 1, 2: 1, 4: 6, 3: 8, 6: 8}
    octa = {**tet, 4: 18, 8: 12}
    oct_sign = {**tet, 4: -6, 8: -12}
    ico = {1: 1, 2: 1, 4: 30, 3: 20, 6: 20, 5: 24, 10: 24}
    groups = [
        ("2T", 12, 24, tet), ("2O", 24, 48, octa),
        ("2O_sign", 24, 48, oct_sign), ("2I", 60, 120, ico),
        ("Q8", 4, 8, {1: 1, 2: 1, 4: 6}),
    ]
    groups.extend(("Q8_" + axis, 4, 8, {1: 1, 2: 1, 4: -2}) for axis in "ijk")
    for name, period, order, signed_counts in groups:
        yield name, period, partial(binary_dimension, order, signed_counts)


def project(period, dim):
    plus, minus, eta, _ = quasipolynomial(period, dim)
    return plus, minus, eta


def shell_record(name, period, dim):
    """One row of data/reference/shells.json: signed multiplicity quasipolynomial and eta."""
    plus, minus, eta, polynomials = quasipolynomial(period, dim)
    degree = max((k for row in polynomials for k, x in enumerate(row) if x), default=0)
    return dict(name=name, period=period, degree=degree,
                ascending_polynomials=[[str(c) for c in row] for row in polynomials],
                eta=str(eta), m_plus=plus, m_minus=minus)


def quasipolynomial(period, dim):
    plus, minus = dim(0, 1), dim(1, 0)
    delta = [dim(n, n + 1) - dim(n + 1, n) for n in range(4 * period)]
    polynomials = []
    for r in range(period):
        v0, v1, v2, v3 = delta[r::period]
        quad = Fraction(v2 - 2*v1 + v0, 2*period**2)
        linear = Fraction(v1 - v0, period) - quad*(2*r + period)
        const = v0 - linear*r - quad*r**2
        if const + linear*(r+3*period) + quad*(r+3*period)**2 != v3:
            raise ArithmeticError("Failure of the angular quasipolynomial")
        polynomials.append((const, linear, quad))
    if any(sum(column) for column in zip(*polynomials)):
        raise ArithmeticError("Non-vanishing average angular coefficient")
    numerator = [
        sum((-1)**j * comb(3, j) * delta[i-j*period] for j in range(i//period+1))
        for i in range(3*period)
    ]
    moments = [sum(comb(i, k)*v for i, v in enumerate(numerator) if i >= k) for k in range(4)]
    if any(moments[:3]):
        raise ArithmeticError("Uncancelled pole of the angular generating function")
    return plus, minus, Fraction(-moments[3], period**3), polynomials
