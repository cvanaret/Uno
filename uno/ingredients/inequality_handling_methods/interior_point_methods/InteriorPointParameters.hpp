// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INTERIORPOINTPARAMETERS_H
#define UNO_INTERIORPOINTPARAMETERS_H

struct InteriorPointParameters {
   double tau_min;
   double k_sigma;
   double dual_regularization_exponent;
   double small_direction_factor;
   double push_variable_to_interior_k1;
   double push_variable_to_interior_k2;
   double damping_factor; // (Section 3.7 in IPOPT paper)
   double default_multiplier;
};

#endif // UNO_INTERIORPOINTPARAMETERS_H