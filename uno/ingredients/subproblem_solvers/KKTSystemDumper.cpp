// Copyright (c) 2026 Charlie Vanaret and Alexis Montoison
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "KKTSystemDumper.hpp"
#include <algorithm>
#include <cctype>
#include <cerrno>
#include <fstream>
#include <sys/stat.h>
#ifdef _WIN32
#include <direct.h>
#endif
#include "COOLinearSystem.hpp"
#include "options/Options.hpp"
#include "symbolic/Range.hpp"
#include "tools/Logger.hpp"

namespace uno {
   KKTSystemDumper::KKTSystemDumper(const Options& options):
         path(options.get_string("dump_kkt_path")),
         // a frequency of 0 would divide by zero; clamp it to "every iteration"
         frequency(std::max<size_t>(1, options.get_unsigned_int("dump_kkt_frequency"))) {
   }

   namespace {
      // replace filesystem-unfriendly characters (e.g. path separators in a model name) with '_'
      std::string sanitize(const std::string& name) {
         std::string result = name;
         for (char& character: result) {
            if (!std::isalnum(static_cast<unsigned char>(character)) && character != '-' && character != '.') {
               character = '_';
            }
         }
         return result;
      }

      // filesystem-friendly base name for the model: drop Uno's reformulation suffix
      // (e.g. " -> bounds relaxed") and any file extension (e.g. ".nl"), then sanitize
      std::string base_name(const std::string& name) {
         std::string base = name.substr(0, name.find(" -> "));
         const size_t dot = base.find_last_of('.');
         if (dot != std::string::npos) {
            base.resize(dot);
         }
         return sanitize(base);
      }

      // create a single directory; success if it now exists (created or already there).
      // std::filesystem is deliberately avoided: it is unavailable on the old macOS
      // deployment targets used by the binary builds
      bool make_directory(const std::string& directory) {
#ifdef _WIN32
         const int return_code = _mkdir(directory.c_str());
#else
         const int return_code = ::mkdir(directory.c_str(), 0777);
#endif
         return return_code == 0 || errno == EEXIST;
      }

      // create a directory and all of its parents (mkdir -p)
      bool create_directories(const std::string& path) {
         for (size_t index = 1; index < path.size(); ++index) {
            if (path[index] == '/' || path[index] == '\\') {
               if (!make_directory(path.substr(0, index))) {
                  return false;
               }
            }
         }
         return make_directory(path);
      }

      // MatrixMarket coordinate/symmetric writer for the lower-triangular COO augmented matrix
      void write_matrix(const COOLinearSystem& system, const std::string& path) {
         std::ofstream stream(path);
         stream << "%%MatrixMarket matrix coordinate real symmetric\n";
         stream << "% regularized KKT (augmented) matrix, lower triangle; duplicate diagonal entries sum\n";
         stream << system.dimension << ' ' << system.dimension << ' ' << system.number_nonzeros << '\n';
         stream.precision(16);
         for (size_t nonzero_index: Range(system.number_nonzeros)) {
            // convert from the solver's index base to MatrixMarket's 1-based indexing
            const auto row = system.matrix_row_indices[nonzero_index] - system.solver_indexing + 1;
            const auto column = system.matrix_column_indices[nonzero_index] - system.solver_indexing + 1;
            stream << row << ' ' << column << ' ' << system.matrix_values[nonzero_index] << '\n';
         }
      }

      // MatrixMarket array writer for a dense right-hand side vector
      void write_vector(const Vector<double>& vector, size_t dimension, const std::string& path) {
         std::ofstream stream(path);
         stream << "%%MatrixMarket matrix array real general\n";
         stream << dimension << " 1\n";
         stream.precision(16);
         for (size_t index: Range(dimension)) {
            stream << vector[index] << '\n';
         }
      }
   }

   void KKTSystemDumper::dump(const LinearSystem& system, const std::string& instance, bool feasibility_phase,
         size_t iteration) const {
      // Only the optimality-phase systems are dumped for now: feasibility-restoration systems are
      // intermittent and do not fit the major-iteration dump frequency (dedicated _feas/_lsq/_soc
      // options may be added later). An empty path disables the dump.
      // The frequency samples the first optimality iteration and every `frequency`-th one after it
      // (iteration is >= 1 here, so the subtraction does not underflow).
      if (feasibility_phase || this->path.empty() || (iteration - 1) % this->frequency != 0) {
         return;
      }
      // the augmented matrix is only available in COO form
      const auto* coo_system = dynamic_cast<const COOLinearSystem*>(&system);
      if (coo_system == nullptr) {
         return;
      }

      if (!create_directories(this->path)) {
         WARNING << "KKTSystemDumper: could not create directory " << this->path << '\n';
         return;
      }

      const std::string stem = this->path + "/kkt_" + base_name(instance) + "_it" + std::to_string(iteration);
      write_matrix(*coo_system, stem + ".mtx");
      write_vector(coo_system->rhs, coo_system->dimension, stem + "_rhs.mtx");
      DEBUG << "KKTSystemDumper: wrote " << stem << ".mtx\n";
   }
} // namespace
