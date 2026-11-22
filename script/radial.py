"""Independent ODE checks of the separated boundary Wronskians."""
from math import acosh, cosh, hypot, isclose, log, sqrt, tanh

from scipy.integrate import solve_ivp


def shell_wronskian(lam, sign, h, radius):
    """Shoot A* A y = 0, with a unit Dirichlet or Robin normalization."""
    t_end = acosh(radius**2) / 2
    xstar = 1.5 / hypot(1, h)
    init = (0., 1.) if lam > 0 else (1., -(lam-sign*xstar))

    def derivative(t, state):
        q = sqrt(cosh(2*t))
        shift = xstar / q**2
        b = (lam-sign*shift)/q
        db_dt = (3*sign*shift-lam)*tanh(2*t)/q
        y, dy_dtau = state
        return q*dy_dtau, (q*b*b-db_dt)*y

    result = solve_ivp(derivative, (0, t_end), init, method="DOP853",
                       rtol=3e-12, atol=1e-14)
    if not result.success:
        raise RuntimeError(result.message)
    y_end, dy_end = result.y[:, -1]
    b_end = (lam-sign*xstar/radius**2)/radius
    return dy_end + b_end*y_end if lam > 0 else y_end


def compare_wronskians(h, radius, J):
    for lam in (-2.5, -1.5, 1.5, 2.5):
        plus = shell_wronskian(lam, 1, h, radius)
        minus = shell_wronskian(lam, -1, h, radius)
        value = log(plus/minus)
        expected = -2*J if lam > 0 else 2*J
        if not isclose(value, expected, abs_tol=3e-10):
            raise ArithmeticError(f"Wronskian discrepancy: {h}, {radius}, {lam}")
