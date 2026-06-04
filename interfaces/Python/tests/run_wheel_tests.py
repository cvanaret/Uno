import unopy

unopy.test_libkrylov()   # raises on failure

# run the example end-to-end
import runpy, os
here = os.path.dirname(__file__)
runpy.run_path(os.path.join(here, "..", "example", "example_hs015.py"), run_name="__main__")