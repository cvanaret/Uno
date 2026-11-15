# Copyright (c) 2026 Alexis Montoison and Charlie Vanaret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

# libkrylov (C interface of Krylov.jl) is loaded at runtime. The test is skipped if no libkrylov is available
# (no bundled library, no UNO_KRYLOV_LIBRARY), unless UNOPY_REQUIRE_LIBKRYLOV=1 (wheel CI).

import os
import sys
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
	assert solver.is_libkrylov_functional(), unopy.libkrylov_error()


def test_libkrylov_path_option_type(solver):
	with pytest.raises(TypeError):
		solver.set_option("libkrylov_path", 1)
	solver.set_option("libkrylov_path", "")


@pytest.mark.skipif(os.name == "nt", reason="POSIX file status flags")
def test_libkrylov_leaves_standard_streams_blocking():
	# the Julia runtime sets O_NONBLOCK on the standard streams, which are shared with the parent process (e.g. the
	# shell of cibuildwheel's container, which exits when its stdin becomes nonblocking). Uno restores the flags
	if not libkrylov_available() and os.environ.get("UNOPY_REQUIRE_LIBKRYLOV") != "1":
		pytest.skip("libkrylov is not available")
	import fcntl
	import subprocess
	read_end, write_end = os.pipe()
	try:
		code = "import unopy; assert unopy.UnoSolver().is_libkrylov_functional()"
		completed = subprocess.run([sys.executable, "-c", code], stdin=read_end, stdout=write_end, check=False)
		assert completed.returncode == 0
		assert not fcntl.fcntl(read_end, fcntl.F_GETFL) & os.O_NONBLOCK
		assert not fcntl.fcntl(write_end, fcntl.F_GETFL) & os.O_NONBLOCK
	finally:
		os.close(read_end)
		os.close(write_end)
