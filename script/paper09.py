"""Recompute tables and render the manuscript figure from the public inputs."""
import argparse
import csv
import json
import math
from pathlib import Path

from scipy.integrate import quad

from angular import lifts, project, shell_record
from radial import compare_wronskians
from robin import signs, record

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "data" / "reference"


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


# Critical Robin brackets of Appendix B: (n, parity) -> (lower, upper) outer radius R/q*.
BRACKETS = {
    "2": {"even": ("1.413463", "1.414964"), "odd": ("1.268809", "1.270310")},
    "3": {"even": ("1.195534", "1.197035"), "odd": ("1.124728", "1.126229")},
    "4": {"even": ("1.116586", "1.118087"), "odd": ("1.073819", "1.075320")},
    "12": {"even": ("1.015053", "1.016554"), "odd": ("1.009276", "1.010674")},
}


def robin_certificate():
    """Directed enclosures at both ends of every critical bracket and at the two sample radii."""
    crit = {n: {name: {"endpoints": {r: record(int(n), r, (parity,)) for r in pair}}
                for parity, (name, pair) in enumerate(cases.items())} for n, cases in BRACKETS.items()}
    radii = {"1.1": {str(n): record(n, "1.1") for n in (2, 3, 4, 5)}, "1.5": {"2": record(2, "1.5")}}
    return {"precision_bits": 112,
            "quantity": "enclosure of m-b at the outer end, a=1, proper-distance L2(dt) normalization; index 1 means m-b<0",
            "critical": crit, "radii": radii}


def recompute(out):
    out.mkdir(parents=True, exist_ok=True)
    shells = [shell_record(name, period, dim) for name, period, dim in lifts()]
    if shells != json.loads((REFERENCE/"shells.json").read_text(encoding="utf-8")):
        raise ValueError("Angular quasipolynomials differ from data/reference/shells.json")
    (out / "shells.json").write_text(json.dumps(shells, indent=1) + "\n", encoding="utf-8")
    cert = robin_certificate()
    if cert != json.loads((REFERENCE/"robin_certificate.json").read_text(encoding="utf-8")):
        raise ValueError("Robin enclosures differ from data/reference/robin_certificate.json")
    (out / "robin_certificate.json").write_text(json.dumps(cert, indent=1) + "\n", encoding="utf-8")
    reference = {r["name"]: r for r in shells}
    rows, ratios, windows, morse = [], [], [], []
    checked = {}
    for n, cases in cert["critical"].items():
        for name, block in cases.items():
            pair = sorted(block["endpoints"], key=float)
            got = [block["endpoints"][radius]["parities"][name]["index"] for radius in pair]
            if got != [1, 0]:
                raise ValueError(f"Failed critical bracket: {n}, {name}, {pair}")
    for radius, degrees in (("1.1", (2, 3, 4, 5)), ("1.5", (2,))):
        checked[radius] = {n: signs(n, radius) for n in degrees}
    for name, period, dim in lifts():
        mplus, mminus, eta = project(period, dim)
        ref = reference[name]
        if (mplus, mminus, str(eta)) != (ref["m_plus"], ref["m_minus"], ref["eta"]):
            raise ValueError(f"Angular mismatch: {name}")
        rows.append(dict(lift=name, m_plus=mplus, m_minus=mminus, eta=str(eta)))
        for h in (0.1, 0.5, 1, 2):
            xstar = 1.5 / math.hypot(1, h)
            u0 = math.asin(xstar / 1.5)/2
            u1 = math.asin(xstar / 2.5)/2
            for window, prep in (("W0", (u0+math.pi/4)/2), ("W1", (u1+u0)/2)):
                for s in (1, -1):
                    shift = xstar/math.sin(2*prep)
                    # k>=1 cannot cross inside these two windows.
                    diff = sum(m*(int(lam-s*shift < 0)-int(lam < 0))
                               for m, lam in ((mplus, 1.5), (mminus, -1.5)))
                    expected = 0 if window == "W0" else s*(mplus if s > 0 else mminus)
                    if diff != expected:
                        raise ArithmeticError("Occupation comparison")
                    windows.append(dict(lift=name, h=h, window=window, s=s, Weyl=diff,
                                        Dirac_vector=0, Dirac_axial=2*diff))
            for rad in (2.0, 3.0, math.inf):
                cutoff = math.acosh(rad**2)/2
                J = quad(lambda t: 2 * xstar*math.exp(-2*t)/(1+math.exp(-4*t)), 0, cutoff)[0]
                closed = 0.5*xstar*math.acos(rad**-2)
                if not math.isclose(J, closed, rel_tol=1e-12):
                    raise ArithmeticError("Radial integral differs from closed form")
                if name == "S3" and math.isfinite(rad):
                    compare_wronskians(h, rad, J)
                delta = -2 * float(eta) * J
                ratios.append(dict(lift=name, h=h, R_over_qstar=rad, J=J, Delta=delta,
                                   modulus_ratio=math.exp(delta), half_weight=math.exp(delta/2)))
        if name in ("S3", "C3_left", "RP3_even", "L5", "L7_23", "2I"):
            for rad, by_n in checked.items():
                ev = 1 + sum(dim(n, n)*idx[0] for n, idx in by_n.items())
                od = sum(dim(n, n)*idx[1] for n, idx in by_n.items())
                morse.append(dict(space=name, R_over_qstar=rad, E1_even=0, E1_odd=0, E2_even=ev, E2_odd=od))
    contacts = []
    for rad in (2, 3):
        T = math.acosh(rad**2)/2
        def density(t):
            B = 1.5*t - 0.75*math.atan(math.sinh(2*t))
            return math.sqrt(math.cosh(2*t))*math.exp(-2*B)
        norm = quad(density, 0, T, epsabs=1e-13)[0]
        contacts.append(dict(R_over_qstar=rad, mu_I=norm, log_mu_I=math.log(norm)))
    out.mkdir(parents=True, exist_ok=True)
    for filename, records in (("angular",rows), ("ratios",ratios), ("windows",windows), ("morse",morse), ("contact",contacts)):
        write_csv(out / (filename+".csv"), records)
    print(f"Verified {len(rows)} lifts, {len(ratios)} relative moduli, {len(windows)} window comparisons, 32 radial Wronskian ratios and eight critical Robin brackets.")


def render(out):
    import matplotlib as mpl
    mpl.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    out.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.25), constrained_layout=True)
    blue, orange, green = "#1d4766", "#b14827", "#267461"
    for data, sign, func, label, col in ((np.linspace(0,.67,220),-1,np.cos,"Lorentz",blue),
                                         (np.linspace(0,.67,220),1,np.cosh,"Euclidean",orange)):
        axes[0].plot(sign*data, func(2*data), label=label, color=col, lw=2.4)
    axes[0].axvline(0, color="0.4", ls=":", lw=1)
    axes[0].set(xlabel="Conformal distance from throat", ylabel=r"$q^2/q_*^2$", xlim=(-.68,.68))
    axes[0].legend(frameon=False)
    u = np.linspace(.18,math.pi/4,301)
    for h, col, style in ((0,blue,"-"),(.5,orange,"--"),(1,green,"-.")):
        axes[1].plot(u,1.5/np.hypot(1,h)/np.sin(2*u), label=f"h = {h}", color=col, ls=style, lw=2)
    for level in (1.5,2.5):
        axes[1].axhline(level, ls=":", lw=1, color="0.35")
    axes[1].set(xlabel="Lorentzian conformal time $u$", ylabel="Axial shift $X_L$", xlim=(.18,math.pi/4+.006), ylim=(.95,3.7))
    axes[1].legend(frameon=False)
    for ax in axes:
        ax.spines[['top','right']].set_visible(False)
    fig.savefig(out/"throat.pdf", metadata={"Creator": "Paper 09 reproducibility package", "CreationDate": None, "ModDate": None})
    fig.savefig(out/"throat.png", dpi=240)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("recompute", "render"))
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    (recompute if args.command == "recompute" else render)(args.output)
