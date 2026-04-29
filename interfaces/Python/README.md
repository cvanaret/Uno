## Uno

Uno is a solver for nonlinearly constrained optimization that unifies Lagrange-Newton (essentially **SQP** and **interior-point**) methods.
It breaks them down into a set of common building blocks (e.g., strategies to compute descent directions and techniques to enforce globalization).

`unopy`, Uno's Python interface, allows you to solve an optimization model described by callback functions.

## Example

An implementation example of the Hock-Schittkowski model [hs015](https://vanderbei.princeton.edu/ampl/nlmodels/hs/hs015.mod) is available in the file [example_hs015.py](https://github.com/cvanaret/Uno/blob/main/interfaces/Python/example/example_hs015.py).

### Querying the current Uno version

Query the current Uno version with:
```python
unopy.current_uno_version()
```

## How to cite Uno

Our Uno paper was accepted in the Mathematical Programming Computation journal on Feb 22, 2026.
Our preprint is available on [ResearchGate](https://www.researchgate.net/publication/397446552_Implementing_a_unified_solver_for_nonlinearly_constrained_optimization) and can be cited with the following BibTeX entry:

```
@unpublished{VanaretLeyffer2026,
  author = {Vanaret, Charlie and Leyffer, Sven},
  title = {Implementing a unified solver for nonlinearly constrained optimization},
  year = {2026},
  note = {Accepted to Mathematical Programming Computation on Feb 22, 2026}
}
```