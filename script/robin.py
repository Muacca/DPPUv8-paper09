"""Validated signs of zero-energy Dirichlet-to-Neumann coefficients."""
from flint import arb, ctx


def riccati_step(value, potential, invp, dt):
    scale = (potential / invp).sqrt()
    th = (dt * (potential * invp).sqrt()).tanh()
    if value is None:
        return scale / th
    return (value + scale * th) / (1 + value * th / scale)


def enclose(radius, n, count):
    r = arb(radius)
    dt = (r**2).acosh() / (2*count)
    k = n*(n+2)
    j = k-3
    even, odd = [arb(0), arb(0)], [None, None]
    left = arb(1)
    for i in range(count):
        right = (2*(i+1)*dt).cosh()
        x = left.union(right)
        invp = k*(18/x+2*j*x)/j
        pot = x*j*(2*j*x**2-6)/(18+2*j*x**2)**2
        for pair in (even, odd):
            pair[0] = riccati_step(pair[0], pot.lower(), invp.upper(), dt).lower()
            pair[1] = riccati_step(pair[1], pot.upper(), invp.lower(), dt).upper()
        left = right
    b = 3 / ((18 + 2*j*r**4)*(r**4-1).sqrt())
    return [(pair[0]-b, pair[1]-b) for pair in (even, odd)]


def resolve(n, radius, parities=(0, 1)):
    """Enclosures of m-b for the requested parities, refined until every sign is decided."""
    ctx.prec = 112
    for count in (128, 512, 2048, 8192, 32768):
        bounds = enclose(radius, n, count)
        if all(bounds[i][1] < 0 or bounds[i][0] > 0 for i in parities):
            return count, bounds
    raise ArithmeticError(f"Unresolved Robin sign for n={n}, R/q*={radius}")


def signs(n, radius, parities=(0, 1)):
    """1 where m-b<0 (a negative direction), 0 where m-b>0."""
    _, bounds = resolve(n, radius, parities)
    return [int(bounds[i][1] < 0) for i in parities]


def record(n, radius, parities=(0, 1)):
    """Serializable enclosure record for data/reference/robin_certificate.json."""
    count, bounds = resolve(n, radius, parities)
    names = ("even", "odd")
    return {"steps": count, "parities": {names[i]: {"lower": bounds[i][0].str(30), "upper": bounds[i][1].str(30),
                                                   "index": int(bounds[i][1] < 0)} for i in parities}}
