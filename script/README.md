# Reproduce Paper 09

Run these commands from the repository root with Python 3.11 or later.

~~~sh
python -m pip install -r script/requirements.txt
python -B script/paper09.py recompute --output results_new
python -B script/paper09.py render --output figures_new
~~~

The first command after installation produces five tables. It starts with signed group-element counts and cyclic weights, recomputes the fifteen eta invariants and lowest-shell multiplicities, and rebuilds the angular reference data (shells.json). A separate directed-ball Riccati integration encloses the Robin quantity at the two endpoints of each of eight critical brackets and at the sample radii needed for the total finite-radius indices (robin_certificate.json). Both reference files are written to the output directory and compared with the distributed copies. A direct second-order ODE shooting check confirms 32 separated Wronskian log ratios at both charge signs. It also computes the first two vacuum windows and relative moduli (including the prescribed half-line limit). Rendering produces PDF and PNG versions of the Lorentzian to Euclidean matching and level crossings from the analytical background. The script writes its named outputs to the chosen output directory. Use a separate directory for comparison against the distributed data.

Tested with Python 3.12.10, NumPy 2.3.5, SciPy 1.16.3, python-flint 0.8.0 and Matplotlib 3.10.7. See the [data inventory](../data/README.md) and [scope of numerical checks](../data/SCOPE.md).
