// Copyright (c) 2026 Charlie Vanaret and Alexis Montoison
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "KKTSystemDumper.hpp"
#include <cctype>
#include <cstdio>
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
         path(options.get_string_optional("dump_kkt_path").value_or("")),
         frequency(options.get_unsigned_int("dump_kkt_frequency")) {
      // a frequency of 0 would divide by zero
      if (this->frequency == 0) {
         throw std::runtime_error("dump_kkt_frequency should be positive");
      }
   }

   NoKKTSystemDumper::NoKKTSystemDumper(const Options& options) {
      if (!options.get_string_optional("dump_kkt_path").value_or("").empty()) {
         WARNING << "'dump_kkt_path' is set, but KKT system dumping is only supported without low-rank corrections\n";
      }
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
      // (e.g. " -> bounds relaxed"), the directory (the AMPL model name is the path of the .nl file)
      // and any file extension (e.g. ".nl"), then sanitize
      std::string base_name(const std::string& name) {
         std::string base = name.substr(0, name.find(" -> "));
         const size_t separator = base.find_last_of("/\\");
         if (separator != std::string::npos) {
            base.erase(0, separator + 1);
         }
         // search the extension after the directory was removed (directories may contain dots)
         const size_t dot = base.find_last_of('.');
         if (dot != std::string::npos) {
            base.resize(dot);
         }
         return sanitize(base);
      }

      bool is_separator(char character) {
         return character == '/' || character == '\\';
      }

      // whether the path exists and is a directory
      bool is_directory(const std::string& directory) {
#ifdef _WIN32
         struct _stat info{};
         return ::_stat(directory.c_str(), &info) == 0 && (info.st_mode & _S_IFDIR) != 0;
#else
         struct stat info{};
         return ::stat(directory.c_str(), &info) == 0 && S_ISDIR(info.st_mode);
#endif
      }

      // create a single directory; success if it now exists (created or already there).
      // std::filesystem is deliberately avoided: it is unavailable on the old macOS
      // deployment targets used by the binary builds
      bool make_directory(const std::string& directory) {
         // check existence first: depending on the platform, mkdir on an existing directory may report
         // EACCES/EROFS (non-writable parent, read-only mount, Windows drive root) instead of EEXIST
         if (is_directory(directory)) {
            return true;
         }
#ifdef _WIN32
         const int return_code = ::_mkdir(directory.c_str());
#else
         const int return_code = ::mkdir(directory.c_str(), 0777);
#endif
         // on failure, the directory may have been created concurrently: check again instead of trusting errno
         return return_code == 0 || is_directory(directory);
      }

      // index of the first character after the root of the path: leading separators (POSIX "/"),
      // a Windows drive ("C:", "C:\\") or a UNC prefix ("\\\\server\\share\\"). The root itself is never created
      size_t root_length(const std::string& path) {
         size_t index = 0;
#ifdef _WIN32
         if (path.size() >= 2 && is_separator(path[0]) && is_separator(path[1])) {
            // UNC path: skip "\\\\server\\share"
            index = 2;
            for (size_t component = 0; component < 2; component) {
               while (index < path.size() && !is_separator(path[index])) {
                  index;
               }
               if (component == 0) {
                  while (index < path.size() && is_separator(path[index])) {
                     index;
                  }
               }
            }
         }
         else if (path.size() >= 2 && std::isalpha(static_cast<unsigned char>(path[0])) && path[1] == ':') {
            index = 2;
         }
#endif
         while (index < path.size() && is_separator(path[index])) {
            ++index;
         }
         return index;
      }

      // create a directory and all of its parents (mkdir -p)
      bool create_directories(const std::string& path) {
         for (size_t index = root_length(path); index < path.size(); ++index) {
            // create each intermediate component once (skip repeated separators, e.g. "a//b")
            if (is_separator(path[index]) && !is_separator(path[index - 1])) {
               if (!make_directory(path.substr(0, index))) {
                  return false;
               }
            }
         }
         return make_directory(path);
      }

      // MatrixMarket coordinate/symmetric writer for the lower-triangular COO augmented matrix
      // returns false if the file could not be opened or written (e.g. permission denied, disk full)
      bool write_matrix(const COOLinearSystem& system, const std::string& path) {
         std::ofstream stream(path);
         if (!stream) {
            return false;
         }
         stream << "%%MatrixMarket matrix coordinate real symmetric\n";
         stream << "% regularized KKT (augmented) matrix, lower triangle; duplicate diagonal entries sum\n";
         stream << system.dimension << ' ' << system.dimension << ' ' << system.number_nonzeros << '\n';
         stream.precision(16);
         for (size_t nonzero_index: Range(system.number_nonzeros)) {
            // convert from the solver's index base to MatrixMarket's 1-based indexing
            const auto row = system.matrix_row_indices[nonzero_index] - system.get_solver_indexing() + 1;
            const auto column = system.matrix_column_indices[nonzero_index] - system.get_solver_indexing() + 1;
            stream << row << ' ' << column << ' ' << system.matrix_values[nonzero_index] << '\n';
         }
         // close explicitly: buffered data is only flushed (and write errors detected) here
         stream.close();
         return !stream.fail();
      }

      // MatrixMarket array writer for a dense right-hand side vector
      // returns false if the file could not be opened or written
      bool write_vector(const Vector<double>& vector, size_t dimension, const std::string& path) {
         std::ofstream stream(path);
         if (!stream) {
            return false;
         }
         stream << "%%MatrixMarket matrix array real general\n";
         stream << dimension << " 1\n";
         stream.precision(16);
         for (size_t index: Range(dimension)) {
            stream << vector[index] << '\n';
         }
         stream.close();
         return !stream.fail();
      }
   }

   void KKTSystemDumper::dump(const LinearSystem& system, const std::string& instance, bool feasibility_phase,
         size_t iteration) const {
      // Only the optimality-phase systems are dumped for now: feasibility-restoration systems are
      // intermittent and do not fit the major-iteration dump frequency (dedicated _feas/_lsq/_soc
      // options may be added later). An empty path disables the dump.
      // The frequency samples the first optimality iteration and every `frequency`-th one after it
      // (iteration is >= 1 here, so the subtraction does not underflow).
      if (this->path.empty() || feasibility_phase || (iteration - 1) % this->frequency != 0) {
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
      const std::string matrix_file = stem + ".mtx";
      const std::string rhs_file = stem + "_rhs.mtx";
      // remove partial output so that a truncated file is never mistaken for a valid dump
      if (!write_matrix(*coo_system, matrix_file)) {
         WARNING << "KKTSystemDumper: could not write " << matrix_file << '\n';
         std::remove(matrix_file.c_str());
         return;
      }
      if (!write_vector(coo_system->rhs, coo_system->dimension, rhs_file)) {
         WARNING << "KKTSystemDumper: could not write " << rhs_file << '\n';
         std::remove(matrix_file.c_str());
         std::remove(rhs_file.c_str());
         return;
      }
      DEBUG << "KKTSystemDumper: wrote " << matrix_file << " and " << rhs_file << '\n';
   }
} // namespace