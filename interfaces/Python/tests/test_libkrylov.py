# Copyright (c) 2026 Alexis Montoison and Charlie Vanaret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

# libkrylov (C interface of Krylov.jl) is loaded at runtime. The test is skipped if no libkrylov is available
# (no bundled library, no UNO_KRYLOV_LIBRARY), unless UNOPY_REQUIRE_LIBKRYLOV=1 (wheel CI).

import os
import pytest
import unopy


def libkrylov_available():
	return unopy.bundled_libkrylov_path() is not None or bool(os.environ.get("UNO_KRYLOV_LIBRARY"))


def test_bundled_libkrylov_present():
	if os.environ.get("UNOPY_REQUIRE_LIBKRYLOV") != "1":
		pytest.skip("UNOPY_REQUIRE_LIBKRYLOV is not set")
	assert unopy.bundled_libkrylov_path() is not None


def test_libkrylov_functional(solver):
	if not libkrylov_available() and os.environ.get("UNOPY_REQUIRE_LIBKRYLOV") != "1":
		pytest.skip("libkrylov is not available")
	assert solver.is_libkrylov_functional()


def test_libkrylov_path_option_type(solver):
	with pytest.raises(TypeError):
		solver.set_option("libkrylov_path", 1)
	solver.set_option("libkrylov_path", "")
