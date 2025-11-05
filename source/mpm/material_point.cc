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

#include <deal.II/mpm/material_point.h>

DEAL_II_NAMESPACE_OPEN

namespace MPM
{
  template <int dim, int spacedim>
  MaterialPoint<dim, spacedim>::MaterialPoint()
    : Particles::Particle<dim, spacedim>()
    , mass(0.0)
    , volume(0.0)
    , velocity()
    , stress()
    , strain()
  {}



  template <int dim, int spacedim>
  MaterialPoint<dim, spacedim>::MaterialPoint(
    const Point<spacedim>      &location,
    const Point<dim>           &reference_location,
    const types::particle_index id)
    : Particles::Particle<dim, spacedim>(location, reference_location, id)
    , mass(0.0)
    , volume(0.0)
    , velocity()
    , stress()
    , strain()
  {}



  template <int dim, int spacedim>
  double
  MaterialPoint<dim, spacedim>::get_mass() const
  {
    return mass;
  }



  template <int dim, int spacedim>
  void
  MaterialPoint<dim, spacedim>::set_mass(const double new_mass)
  {
    mass = new_mass;
  }



  template <int dim, int spacedim>
  double
  MaterialPoint<dim, spacedim>::get_volume() const
  {
    return volume;
  }



  template <int dim, int spacedim>
  void
  MaterialPoint<dim, spacedim>::set_volume(const double new_volume)
  {
    volume = new_volume;
  }



  template <int dim, int spacedim>
  Tensor<1, spacedim>
  MaterialPoint<dim, spacedim>::get_velocity() const
  {
    return velocity;
  }



  template <int dim, int spacedim>
  void
  MaterialPoint<dim, spacedim>::set_velocity(
    const Tensor<1, spacedim> &new_velocity)
  {
    velocity = new_velocity;
  }



  template <int dim, int spacedim>
  SymmetricTensor<2, spacedim>
  MaterialPoint<dim, spacedim>::get_stress() const
  {
    return stress;
  }



  template <int dim, int spacedim>
  void
  MaterialPoint<dim, spacedim>::set_stress(
    const SymmetricTensor<2, spacedim> &new_stress)
  {
    stress = new_stress;
  }



  template <int dim, int spacedim>
  SymmetricTensor<2, spacedim>
  MaterialPoint<dim, spacedim>::get_strain() const
  {
    return strain;
  }



  template <int dim, int spacedim>
  void
  MaterialPoint<dim, spacedim>::set_strain(
    const SymmetricTensor<2, spacedim> &new_strain)
  {
    strain = new_strain;
  }

} // namespace MPM

DEAL_II_NAMESPACE_CLOSE
