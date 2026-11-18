# Copyright (c) 2026 Francois Gallard, Jean-Christophe Giret
# Licensed under the MIT license. See LICENSE file in the project directory for details.
"""A ``scipy.optimize.minimize``-like interface to Uno.

This subpackage requires SciPy, e.g. ``pip install unopy[scipy]``.

Examples
--------
>>> from unopy.scipy import minimize
"""

from __future__ import annotations

from unopy.scipy._minimize import minimize

__all__ = ["minimize"]
