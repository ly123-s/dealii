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

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/mpm/mpm_handler.h>

DEAL_II_NAMESPACE_OPEN

namespace MPM
{
  template <int dim, int spacedim>
  MPMHandler<dim, spacedim>::MPMHandler()
    : particle_handler()
    , dof_handler()
    , fe(nullptr)
    , triangulation(nullptr)
    , mapping(nullptr)
  {}



  template <int dim, int spacedim>
  MPMHandler<dim, spacedim>::MPMHandler(
    const Triangulation<dim, spacedim> &tria,
    const Mapping<dim, spacedim>       &map)
    : particle_handler(tria, map)
    , dof_handler(tria)
    , fe(nullptr)
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
    dof_handler.reinit(tria);
  }



  template <int dim, int spacedim>
  void
  MPMHandler<dim, spacedim>::setup_background_dofs(
    const FiniteElement<dim, spacedim> &finite_element)
  {
    fe = &finite_element;
    dof_handler.distribute_dofs(*fe);
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

    // Check if FE is set up
    if (fe == nullptr)
      return;

    // Create temporary momentum vector
    VectorType grid_momentum;
    grid_momentum.reinit(grid_velocity);

    // Set up quadrature rule for shape function evaluation
    const QGauss<dim> quadrature(fe->degree + 1);
    FEValues<dim, spacedim> fe_values(*mapping,
                                       *fe,
                                       quadrature,
                                       update_values | update_gradients |
                                         update_quadrature_points);

    const unsigned int dofs_per_cell = fe->n_dofs_per_cell();
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    // Loop over all particles
    for (const auto &particle : particle_handler)
      {
        // Get particle properties (assuming: mass, volume, vx, vy, [vz])
        const auto properties = particle.get_properties();
        if (properties.size() < 2 + dim)
          continue;

        const double                mass     = properties[0];
        const Tensor<1, spacedim>   velocity = [&]() {
          Tensor<1, spacedim> vel;
          for (unsigned int d = 0; d < dim; ++d)
            vel[d] = properties[2 + d];
          return vel;
        }();

        // Get the cell containing this particle
        const auto cell = particle.get_surrounding_cell();
        fe_values.reinit(cell);
        cell->get_dof_indices(local_dof_indices);

        // Evaluate shape functions at particle location
        const Point<dim> &ref_location = particle.get_reference_location();
        
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            const double shape_value = fe->shape_value(i, ref_location);
            const types::global_dof_index global_index = local_dof_indices[i];

            // Transfer mass
            grid_mass[global_index] += mass * shape_value;

            // Transfer momentum (for each component)
            for (unsigned int d = 0; d < dim; ++d)
              {
                const types::global_dof_index vel_index = 
                  global_index * dim + d;
                grid_momentum[vel_index] += 
                  mass * velocity[d] * shape_value;
              }
          }
      }

    // Compute velocities from momentum and mass
    for (unsigned int i = 0; i < grid_mass.size(); ++i)
      {
        if (grid_mass[i] > 1e-15)
          {
            for (unsigned int d = 0; d < dim; ++d)
              {
                const types::global_dof_index vel_index = i * dim + d;
                grid_velocity[vel_index] = 
                  grid_momentum[vel_index] / grid_mass[i];
              }
          }
      }
  }



  template <int dim, int spacedim>
  template <typename VectorType>
  void
  MPMHandler<dim, spacedim>::compute_grid_forces(VectorType &grid_force) const
  {
    // Reset grid forces
    grid_force = 0.0;

    // Check if FE is set up
    if (fe == nullptr)
      return;

    // Set up quadrature rule
    const QGauss<dim> quadrature(fe->degree + 1);
    FEValues<dim, spacedim> fe_values(*mapping,
                                       *fe,
                                       quadrature,
                                       update_values | update_gradients);

    const unsigned int dofs_per_cell = fe->n_dofs_per_cell();
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    // Loop over all particles
    for (const auto &particle : particle_handler)
      {
        // Get particle properties
        // Assuming properties: mass, volume, velocities, then stress components
        const auto properties = particle.get_properties();
        if (properties.size() < 2 + dim)
          continue;

        const double volume = properties[1];

        // For now, assume stress is stored after velocity components
        // In a full implementation, stress would be in MaterialPoint
        // Here we use a simplified zero stress for demonstration
        SymmetricTensor<2, spacedim> stress;
        
        // Get cell and evaluate shape function gradients
        const auto cell = particle.get_surrounding_cell();
        fe_values.reinit(cell);
        cell->get_dof_indices(local_dof_indices);

        const Point<dim> &ref_location = particle.get_reference_location();

        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            const Tensor<1, dim> shape_grad = 
              fe->shape_grad(i, ref_location);
            const types::global_dof_index global_index = local_dof_indices[i];

            // Compute internal force: f_i = -Volume * stress * grad(N_i)
            const Tensor<1, spacedim> force = -volume * (stress * shape_grad);

            for (unsigned int d = 0; d < dim; ++d)
              {
                const types::global_dof_index force_index = 
                  global_index * dim + d;
                grid_force[force_index] += force[d];
              }
          }
      }
  }



  template <int dim, int spacedim>
  template <typename VectorType>
  void
  MPMHandler<dim, spacedim>::grid_to_particle(const VectorType &grid_velocity,
                                               const double      dt)
  {
    // Check if FE is set up
    if (fe == nullptr)
      return;

    // Set up quadrature rule
    const QGauss<dim> quadrature(fe->degree + 1);
    FEValues<dim, spacedim> fe_values(*mapping,
                                       *fe,
                                       quadrature,
                                       update_values | update_gradients);

    const unsigned int dofs_per_cell = fe->n_dofs_per_cell();
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    // Loop over all particles
    for (auto &particle : particle_handler)
      {
        // Get particle properties
        auto properties = particle.get_properties();
        if (properties.size() < 2 + dim)
          continue;

        // Get cell
        const auto cell = particle.get_surrounding_cell();
        fe_values.reinit(cell);
        cell->get_dof_indices(local_dof_indices);

        const Point<dim> &ref_location = particle.get_reference_location();

        // Interpolate velocity from grid to particle
        Tensor<1, spacedim> new_velocity;
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            const double shape_value = fe->shape_value(i, ref_location);
            const types::global_dof_index global_index = local_dof_indices[i];

            for (unsigned int d = 0; d < dim; ++d)
              {
                const types::global_dof_index vel_index = 
                  global_index * dim + d;
                new_velocity[d] += grid_velocity[vel_index] * shape_value;
              }
          }

        // Update particle velocity in properties
        for (unsigned int d = 0; d < dim; ++d)
          properties[2 + d] = new_velocity[d];

        // Update particle position
        Point<spacedim> new_position = particle.get_location();
        for (unsigned int d = 0; d < dim; ++d)
          new_position[d] += new_velocity[d] * dt;

        particle.set_location(new_position);
        particle.set_properties(properties);
      }

    // Update particle reference locations after position changes
    particle_handler.sort_particles_into_subdomains_and_cells();
  }



  template <int dim, int spacedim>
  void
  MPMHandler<dim, spacedim>::update_stresses(const double dt)
  {
    // Check if FE is set up
    if (fe == nullptr)
      return;

    // Set up quadrature rule
    const QGauss<dim> quadrature(fe->degree + 1);
    FEValues<dim, spacedim> fe_values(*mapping,
                                       *fe,
                                       quadrature,
                                       update_values | update_gradients);

    const unsigned int dofs_per_cell = fe->n_dofs_per_cell();
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    // Loop over all particles
    for (auto &particle : particle_handler)
      {
        // Get particle properties
        auto properties = particle.get_properties();
        if (properties.size() < 2 + dim)
          continue;

        // Get cell
        const auto cell = particle.get_surrounding_cell();
        fe_values.reinit(cell);
        cell->get_dof_indices(local_dof_indices);

        const Point<dim> &ref_location = particle.get_reference_location();

        // Compute velocity gradient at particle location
        Tensor<2, spacedim> velocity_gradient;
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            const Tensor<1, dim> shape_grad = 
              fe->shape_grad(i, ref_location);
            const types::global_dof_index global_index = local_dof_indices[i];

            Tensor<1, spacedim> velocity;
            for (unsigned int d = 0; d < dim; ++d)
              velocity[d] = properties[2 + d];

            for (unsigned int i_dim = 0; i_dim < dim; ++i_dim)
              for (unsigned int j_dim = 0; j_dim < dim; ++j_dim)
                velocity_gradient[i_dim][j_dim] += 
                  velocity[i_dim] * shape_grad[j_dim];
          }

        // Compute strain rate (symmetric part of velocity gradient)
        SymmetricTensor<2, spacedim> strain_rate;
        for (unsigned int i_dim = 0; i_dim < dim; ++i_dim)
          for (unsigned int j_dim = i_dim; j_dim < dim; ++j_dim)
            strain_rate[i_dim][j_dim] = 
              0.5 * (velocity_gradient[i_dim][j_dim] + 
                     velocity_gradient[j_dim][i_dim]);

        // Compute strain increment
        const SymmetricTensor<2, spacedim> strain_increment = 
          strain_rate * dt;

        // For now, we just store the strain increment
        // In a full implementation, this would update stress using
        // a constitutive model (e.g., LinearElasticModel)
        // and store it back in particle properties
        
        particle.set_properties(properties);
      }
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
    if (!properties.empty())
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
