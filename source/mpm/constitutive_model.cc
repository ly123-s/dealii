// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/mpm/constitutive_model.h>

DEAL_II_NAMESPACE_OPEN

namespace MPM
{
  template <int dim>
  LinearElasticModel<dim>::LinearElasticModel(const double E, const double nu)
    : youngs_modulus(E)
    , poisson_ratio(nu)
  {
    // Compute Lamé parameters from Young's modulus and Poisson's ratio
    lambda = (E * nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
    mu     = E / (2.0 * (1.0 + nu));
  }



  template <int dim>
  void
  LinearElasticModel<dim>::update_stress(
    SymmetricTensor<2, dim>       &stress,
    const SymmetricTensor<2, dim> &strain_increment,
    const double                   dt) const
  {
    (void)dt; // Time step not needed for linear elasticity

    // Compute stress increment using linear elastic constitutive law
    // sigma_ij = lambda * trace(epsilon) * delta_ij + 2 * mu * epsilon_ij
    
    const double trace_strain = trace(strain_increment);
    
    SymmetricTensor<2, dim> stress_increment;
    for (unsigned int i = 0; i < dim; ++i)
      {
        for (unsigned int j = i; j < dim; ++j)
          {
            if (i == j)
              stress_increment[i][j] = lambda * trace_strain + 
                                        2.0 * mu * strain_increment[i][j];
            else
              stress_increment[i][j] = 2.0 * mu * strain_increment[i][j];
          }
      }
    
    // Update stress
    stress += stress_increment;
  }

} // namespace MPM

DEAL_II_NAMESPACE_CLOSE
