// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PROGRESSMEASURES_H
#define UNO_PROGRESSMEASURES_H

struct ProgressMeasures {
   double infeasibility{}; // constraint violation
   double optimality{}; // optimality measure (scalable by penalty parameter): objective, Lagrangian, ...
};

#endif // UNO_PROGRESSMEASURES_H
