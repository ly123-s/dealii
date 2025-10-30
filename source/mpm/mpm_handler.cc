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

#include <deal.II/grid/grid_tools.h>

#include <deal.II/mpm/mpm_handler.h>

DEAL_II_NAMESPACE_OPEN

namespace MPM
{
  template <int dim, int spacedim>
  MPMHandler<dim, spacedim>::MPMHandler()
    : particle_handler()
    , triangulation(nullptr)
    , mapping(nullptr)
  {}



  template <int dim, int spacedim>
  MPMHandler<dim, spacedim>::MPMHandler(
    const Triangulation<dim, spacedim> &tria,
    const Mapping<dim, spacedim>       &map)
    : particle_handler(tria, map)
    , triangulation(&tria)
    , mapping(&map)
  {}



  template <int dim, int spacedim>
  void
  MPMHandler<dim, spacedim>::initialize(
    const Triangulation<dim, spacedim> &tria,
    const Mapping<dim, spacedim>       &map,
    const unsigned int                  n_properties)
  {
    triangulation = &tria;
    mapping       = &map;
    particle_handler.initialize(tria, map, n_properties);
  }



  template <int dim, int spacedim>
  template <typename VectorType>
  void
  MPMHandler<dim, spacedim>::particle_to_grid(VectorType &grid_velocity,
                                               VectorType &grid_mass) const
  {
    // Reset grid quantities
    grid_velocity = 0.0;
    grid_mass     = 0.0;

    // This is a simplified implementation placeholder.
    // In a full implementation, this would:
    // 1. Loop over all material points
    // 2. Compute shape functions at material point positions
    // 3. Transfer mass and momentum to grid nodes using shape functions
    // 4. Compute grid velocities from momentum and mass
  }



  template <int dim, int spacedim>
  template <typename VectorType>
  void
  MPMHandler<dim, spacedim>::compute_grid_forces(VectorType &grid_force) const
  {
    // Reset grid forces
    grid_force = 0.0;

    // This is a simplified implementation placeholder.
    // In a full implementation, this would:
    // 1. Loop over all material points
    // 2. Compute shape function gradients at material point positions
    // 3. Compute internal forces from stress and volume
    // 4. Transfer forces to grid nodes
  }



  template <int dim, int spacedim>
  template <typename VectorType>
  void
  MPMHandler<dim, spacedim>::grid_to_particle(const VectorType &grid_velocity,
                                               const double      dt)
  {
    // This is a simplified implementation placeholder.
    // In a full implementation, this would:
    // 1. Loop over all material points
    // 2. Interpolate grid velocities to material point positions
    // 3. Update material point velocities
    // 4. Update material point positions based on velocities and dt
    (void)grid_velocity;
    (void)dt;
  }



  template <int dim, int spacedim>
  void
  MPMHandler<dim, spacedim>::update_stresses(const double dt)
  {
    // This is a simplified implementation placeholder.
    // In a full implementation, this would:
    // 1. Loop over all material points
    // 2. Compute deformation gradient from grid velocity gradients
    // 3. Compute strain increment
    // 4. Update stress using constitutive model
    (void)dt;
  }



  template <int dim, int spacedim>
  types::particle_index
  MPMHandler<dim, spacedim>::n_global_particles() const
  {
    return particle_handler.n_global_particles();
  }



  template <int dim, int spacedim>
  void
  MPMHandler<dim, spacedim>::insert_particle(
    const Point<spacedim>      &position,
    const std::vector<double>  &properties,
    const types::particle_index id)
  {
    // Find the cell containing the particle position
    auto cell = GridTools::find_active_cell_around_point(*mapping,
                                                          *triangulation,
                                                          position);
    
    // Create a new particle
    Particles::Particle<dim, spacedim> new_particle;
    new_particle.set_location(position);
    new_particle.set_reference_location(
      mapping->transform_real_to_unit_cell(cell, position));
    new_particle.set_id(id);
    
    // Set properties if provided
    if (properties.size() > 0)
      new_particle.set_properties(properties);
    
    // Insert the particle
    particle_handler.insert_particle(new_particle, cell);
  }



  template <int dim, int spacedim>
  void
  MPMHandler<dim, spacedim>::update_cached_numbers()
  {
    particle_handler.update_cached_numbers();
  }



  template <int dim, int spacedim>
  void
  MPMHandler<dim, spacedim>::clear_particles()
  {
    particle_handler.clear_particles();
  }

} // namespace MPM

DEAL_II_NAMESPACE_CLOSE
