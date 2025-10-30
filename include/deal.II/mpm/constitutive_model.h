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

#ifndef dealii_mpm_constitutive_model_h
#define dealii_mpm_constitutive_model_h

#include <deal.II/base/config.h>

#include <deal.II/base/tensor.h>

DEAL_II_NAMESPACE_OPEN

namespace MPM
{
  /**
   * Base class for constitutive models in MPM.
   * This class defines the interface for material constitutive laws
   * that relate stress to strain and deformation.
   */
  template <int dim>
  class ConstitutiveModel
  {
  public:
    /**
     * Virtual destructor.
     */
    virtual ~ConstitutiveModel() = default;

    /**
     * Update stress based on strain increment.
     *
     * @param stress Current stress tensor (input/output).
     * @param strain_increment Incremental strain tensor.
     * @param dt Time step size.
     */
    virtual void
    update_stress(SymmetricTensor<2, dim>       &stress,
                  const SymmetricTensor<2, dim> &strain_increment,
                  const double                   dt) const = 0;
  };


  /**
   * Linear elastic constitutive model.
   * This implements Hooke's law for linear elasticity.
   */
  template <int dim>
  class LinearElasticModel : public ConstitutiveModel<dim>
  {
  public:
    /**
     * Constructor.
     *
     * @param youngs_modulus Young's modulus of the material.
     * @param poisson_ratio Poisson's ratio of the material.
     */
    LinearElasticModel(const double youngs_modulus,
                       const double poisson_ratio);

    /**
     * Update stress based on strain increment using linear elasticity.
     */
    virtual void
    update_stress(SymmetricTensor<2, dim>       &stress,
                  const SymmetricTensor<2, dim> &strain_increment,
                  const double                   dt) const override;

  private:
    /**
     * Young's modulus (E).
     */
    double youngs_modulus;

    /**
     * Poisson's ratio (nu).
     */
    double poisson_ratio;

    /**
     * First Lamé parameter (lambda).
     */
    double lambda;

    /**
     * Second Lamé parameter (mu, shear modulus).
     */
    double mu;
  };

} // namespace MPM

DEAL_II_NAMESPACE_CLOSE

#endif
