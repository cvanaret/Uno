#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace uno {
   // individual classes
   void define_SparseVector(py::module& module);
   void define_RectangularMatrix(py::module& module);
   void define_Options(py::module& module);
   //void define_Model(py::module& module);
   //void define_UnoSolver(py::module& module);

   // module definition
   PYBIND11_MODULE(unopy, module) {
      module.doc() = "Python binding to the solver Uno for nonlinearly constrained optimization";
      
      define_SparseVector(module);
      define_RectangularMatrix(module);
      define_Options(module);
      //define_Model(module);
      //define_UnoSolver(module);
   }
} // namespace
