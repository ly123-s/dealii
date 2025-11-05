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

#ifndef dealii_mpm_material_point_h
#define dealii_mpm_material_point_h

#include <deal.II/base/config.h>

#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/particles/particle.h>

DEAL_II_NAMESPACE_OPEN

/**
 * A namespace for all classes and functions related to the Material Point
 * Method (MPM) implementation.
 */
namespace MPM
{
  /**
   * This class represents a material point in the Material Point Method.
   * It extends the basic Particle concept with material-specific properties
   * such as mass, volume, velocity, stress, and strain.
   *
   * Material points are the fundamental computational objects in MPM,
   * carrying both kinematic and constitutive information as they move
   * through the computational domain.
   */
  template <int dim, int spacedim = dim>
  class MaterialPoint : public Particles::Particle<dim, spacedim>
  {
  public:
    /**
     * Default constructor.
     */
    MaterialPoint();

    /**
     * Constructor that initializes a material point at a given position
     * with specified properties.
     *
     * @param location The position of the material point.
     * @param reference_location The reference position in the element.
     * @param id The unique identifier for this material point.
     */
    MaterialPoint(const Point<spacedim> &location,
                  const Point<dim>      &reference_location,
                  const types::particle_index id);

    /**
     * Get the mass of the material point.
     */
    double
    get_mass() const;

    /**
     * Set the mass of the material point.
     */
    void
    set_mass(const double mass);

    /**
     * Get the volume of the material point.
     */
    double
    get_volume() const;

    /**
     * Set the volume of the material point.
     */
    void
    set_volume(const double volume);

    /**
     * Get the velocity of the material point.
     */
    Tensor<1, spacedim>
    get_velocity() const;

    /**
     * Set the velocity of the material point.
     */
    void
    set_velocity(const Tensor<1, spacedim> &velocity);

    /**
     * Get the Cauchy stress tensor at the material point.
     */
    SymmetricTensor<2, spacedim>
    get_stress() const;

    /**
     * Set the Cauchy stress tensor at the material point.
     */
    void
    set_stress(const SymmetricTensor<2, spacedim> &stress);

    /**
     * Get the strain tensor at the material point.
     */
    SymmetricTensor<2, spacedim>
    get_strain() const;

    /**
     * Set the strain tensor at the material point.
     */
    void
    set_strain(const SymmetricTensor<2, spacedim> &strain);

  private:
    /**
     * Mass of the material point.
     */
    double mass;

    /**
     * Volume of the material point.
     */
    double volume;

    /**
     * Velocity vector.
     */
    Tensor<1, spacedim> velocity;

    /**
     * Cauchy stress tensor.
     */
    SymmetricTensor<2, spacedim> stress;

    /**
     * Strain tensor.
     */
    SymmetricTensor<2, spacedim> strain;
  };

} // namespace MPM

DEAL_II_NAMESPACE_CLOSE

#endif
