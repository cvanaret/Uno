// Copyright (c) 2026 Charlie Vanaret and Alexis Montoison
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_KKTSYSTEMDUMPER_H
#define UNO_KKTSYSTEMDUMPER_H

#include <cstddef>
#include <string>

namespace uno {
   class LinearSystem;
   class Options;

   // Dumps the KKT (augmented) linear system to disk in MatrixMarket format for offline analysis.
   // Controlled by the options dump_kkt_path / dump_kkt_frequency. When disabled (empty path),
   // every call is a cheap no-op.
   class KKTSystemDumper {
   public:
      explicit KKTSystemDumper(const Options& options);

      // Dump the regularized augmented matrix and its right-hand side. `instance` names the model,
      // `feasibility_phase` marks feasibility-restoration systems (currently skipped; only optimality
      // systems are dumped), `iteration` is the major iteration (used for the filename and frequency).
      void dump(const LinearSystem& system, const std::string& instance, bool feasibility_phase, size_t iteration) const;

   private:
      // output directory; an empty path disables the dump
      const std::string path;
      const size_t frequency;
   };
} // namespace

#endif // UNO_KKTSYSTEMDUMPER_H
