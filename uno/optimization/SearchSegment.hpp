// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SEARCHSEGMENT_H
#define UNO_SEARCHSEGMENT_H

namespace uno {
   // forward declaration
   class Iterate;

   class SearchSegment {
      // trial(α) = current ⊕ (base_offset + α · step), optionally projected onto [l, u]
      void assemble_trial_iterate(const Iterate& current, double step_length, Iterate& trial) const;
      [[nodiscard]] double predicted_reduction(double step_length) const;
   };
} // namespace

#endif // UNO_SEARCHSEGMENT_H