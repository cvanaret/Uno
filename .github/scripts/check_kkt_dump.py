"""Validate the KKT systems dumped by Uno (dump_kkt_path) in MatrixMarket format.

usage: check_kkt_dump.py DIR INSTANCE EXPECTED_ITERATIONS NUMBER_CONSTRAINTS
- EXPECTED_ITERATIONS: comma-separated list of iterations that must have been dumped (and only those)
- NUMBER_CONSTRAINTS: expected number of negative eigenvalues of the regularized KKT matrix
  (pass -1 to skip the inertia check, e.g. for ill-conditioned late interior-point iterations)
"""
import math
import pathlib
import re
import sys

import numpy as np

directory, instance = pathlib.Path(sys.argv[1]), sys.argv[2]
expected_iterations = sorted(int(i) for i in sys.argv[3].split(","))
number_constraints = int(sys.argv[4])
errors = []


def fail(message):
    errors.append(message)
    print(f"::error::{message}")


def data_lines(path):
    # skip the header and the comment lines
    return [line.split() for line in path.read_text().splitlines()[1:] if line and not line.startswith("%")]


# 1. file names: exactly one (matrix, rhs) pair per expected iteration, named kkt_<instance>_it<k>
pattern = re.compile(rf"kkt_{re.escape(instance)}_it(\d+)(_rhs)?\.mtx")
found = sorted(p.name for p in directory.iterdir())
iterations = sorted({int(m.group(1)) for name in found if (m := pattern.fullmatch(name))})
unexpected = [name for name in found if not pattern.fullmatch(name)]
if unexpected:
    fail(f"unexpected files (expected kkt_{instance}_it<k>[_rhs].mtx): {unexpected[:4]}")
if iterations != expected_iterations:
    fail(f"dumped iterations {iterations}, expected {expected_iterations}")

for iteration in iterations:
    matrix_path = directory / f"kkt_{instance}_it{iteration}.mtx"
    rhs_path = directory / f"kkt_{instance}_it{iteration}_rhs.mtx"
    if not matrix_path.exists() or not rhs_path.exists():
        fail(f"it{iteration}: missing matrix or RHS file")
        continue

    # 2. matrix: header, size line, 1-based lower-triangular entries, finite values
    if matrix_path.read_text().splitlines()[0] != "%%MatrixMarket matrix coordinate real symmetric":
        fail(f"{matrix_path.name}: wrong header")
    lines = data_lines(matrix_path)
    rows, columns, nonzeros = map(int, lines[0])
    entries = lines[1:]
    if rows != columns:
        fail(f"{matrix_path.name}: not square ({rows}x{columns})")
    if len(entries) != nonzeros:
        fail(f"{matrix_path.name}: header announces {nonzeros} nonzeros, file has {len(entries)}")
    matrix = np.zeros((rows, rows))
    for row, column, value in entries:
        row, column, value = int(row), int(column), float(value)
        if not (1 <= column <= row <= rows):
            fail(f"{matrix_path.name}: entry ({row}, {column}) is out of range or not in the lower triangle")
            break
        if not math.isfinite(value):
            fail(f"{matrix_path.name}: non-finite value at ({row}, {column})")
            break
        # duplicate entries sum (regularization is stored as extra diagonal entries)
        matrix[row - 1, column - 1] += value
    matrix = matrix + np.tril(matrix, -1).T

    # 3. RHS: header, dimension consistent with the matrix, finite values
    if rhs_path.read_text().splitlines()[0] != "%%MatrixMarket matrix array real general":
        fail(f"{rhs_path.name}: wrong header")
    lines = data_lines(rhs_path)
    if list(map(int, lines[0])) != [rows, 1]:
        fail(f"{rhs_path.name}: size {lines[0]}, expected [{rows}, 1]")
    rhs = np.array([float(line[0]) for line in lines[1:]])
    if rhs.size != rows or not np.all(np.isfinite(rhs)):
        fail(f"{rhs_path.name}: {rhs.size} entries (expected {rows}) or non-finite values")
        continue

    # 4. the dumped matrix is the regularized one: it has the correct inertia (n, m, 0) and the system is solvable
    eigenvalues = np.linalg.eigvalsh(matrix)
    tolerance = 1e-12 * max(1., np.abs(eigenvalues).max())
    inertia = (int((eigenvalues > tolerance).sum()), int((eigenvalues < -tolerance).sum()),
               int((np.abs(eigenvalues) <= tolerance).sum()))
    if number_constraints >= 0 and inertia != (rows - number_constraints, number_constraints, 0):
        fail(f"{matrix_path.name}: inertia {inertia}, expected ({rows - number_constraints}, {number_constraints}, 0)")
    solution = np.linalg.solve(matrix, rhs)
    residual = np.linalg.norm(matrix @ solution - rhs) / max(1., np.linalg.norm(rhs))
    if not np.all(np.isfinite(solution)) or residual > 1e-8:
        fail(f"{matrix_path.name}: system not solvable (relative residual {residual:.2e})")

print(f"{instance}: checked {len(iterations)} dumped systems in {directory}")
sys.exit(1 if errors else 0)
