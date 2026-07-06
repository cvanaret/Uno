# Copyright (c) 2026 Charlie Vanaret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

import unopy; print("import OK", flush=True)

# test libkrylov (raises on failure)
unopy.test_libkrylov(); print("call OK", flush=True)

# run the hs015 example
import runpy, os
here = os.path.dirname(__file__)
runpy.run_path(os.path.join(here, "..", "example", "example_hs015.py"), run_name="__main__")