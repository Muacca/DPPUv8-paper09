# Data for Paper 09

| File | Contents |
|---|---|
| [angular.csv](tables/angular.csv) | Exact lowest-shell multiplicities and ordinary Dirac eta invariant for fifteen oriented lifts |
| [ratios.csv](tables/ratios.csv) | 180 fixed-domain log ratios, modulus and half-log weight ratios, for four real nonzero Holst couplings and three outer radii |
| [windows.csv](tables/windows.csv) | Weyl number, Dirac vector and axial differences for both charge signs and two Lorentzian preparation windows |
| [morse.csv](tables/morse.csv) | Full spherical bosonic even/odd index at two representative finite radii, including the homogeneous Fourier direction |
| [contact.csv](tables/contact.csv) | Change-of-domain correction at zero Holst coupling for two finite cutoff radii |
| [shells.json](reference/shells.json) | Signed multiplicity quasipolynomial for each lift, coefficients in ascending powers of the full shell label on each residue class |
| [robin_certificate.json](reference/robin_certificate.json) | Directed enclosures of the Robin quantity at both ends of the eight critical brackets and at the sample radii 1.1 and 1.5, with step counts and signs |

All files here are written by `python -B script/paper09.py recompute`, which also compares each of them with the distributed copy. The two reference files hold the angular projections behind the relative Dirac factor (shells.json) and the directed Robin enclosures behind the real three-form Morse indices (robin_certificate.json), in full precision. Their lifts, scaling and angular conventions are described in the [English PDF](../DPPUv8-paper09.pdf). [Verification scope](SCOPE.md).
