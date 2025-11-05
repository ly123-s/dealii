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

#ifndef dealii_mpm_mpm_handler_h
#define dealii_mpm_mpm_handler_h

#include <deal.II/base/config.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#include <deal.II/mpm/material_point.h>

#include <deal.II/particles/particle_handler.h>

DEAL_II_NAMESPACE_OPEN

namespace MPM
{
  /**
   * This class manages material points for the Material Point Method (MPM).
   * It extends the functionality of the ParticleHandler class with
   * MPM-specific operations such as particle-to-grid and grid-to-particle
   * transfers, stress updates, and constitutive model evaluation.
   *
   * The MPMHandler maintains a collection of MaterialPoint objects and
   * provides methods for:
   * - Transferring momentum and mass from material points to grid nodes
   * - Computing forces on grid nodes from material point stresses
   * - Updating material point positions and velocities from grid solutions
   * - Updating material point stresses based on deformation
   */
  template <int dim, int spacedim = dim>
  class MPMHandler
  {
  public:
    /**
     * Default constructor.
     */
    MPMHandler();

    /**
     * Constructor that initializes the MPM handler with a given
     * triangulation and mapping.
     *
     * @param tria The triangulation representing the background grid.
     * @param mapping The mapping used to transform between reference and
     * real cells.
     */
    MPMHandler(const Triangulation<dim, spacedim> &tria,
               const Mapping<dim, spacedim>       &mapping);

    /**
     * Initialize the MPM handler with a triangulation and mapping.
     *
     * @param tria The triangulation representing the background grid.
     * @param mapping The mapping used to transform between reference and
     * real cells.
     * @param n_properties The number of properties each material point stores.
     */
    void
    initialize(const Triangulation<dim, spacedim> &tria,
               const Mapping<dim, spacedim>       &mapping,
               const unsigned int                  n_properties = 0);

    /**
     * Set up the finite element and DoF handler for grid operations.
     * This is required before using particle_to_grid and grid_to_particle.
     *
     * @param fe The finite element to use for the background grid.
     */
    void
    setup_background_dofs(const FiniteElement<dim, spacedim> &fe);

    /**
     * Transfer mass and momentum from material points to grid nodes.
     * This is the particle-to-grid (P2G) transfer step in MPM.
     *
     * @param grid_velocity Output vector for grid node velocities.
     * @param grid_mass Output vector for grid node masses.
     */
    template <typename VectorType>
    void
    particle_to_grid(VectorType &grid_velocity, VectorType &grid_mass) const;

    /**
     * Compute internal forces on grid nodes from material point stresses.
     * This contributes to the force calculation on the background grid.
     *
     * @param grid_force Output vector for grid node forces.
     */
    template <typename VectorType>
    void
    compute_grid_forces(VectorType &grid_force) const;

    /**
     * Update material point velocities and positions from grid solution.
     * This is the grid-to-particle (G2P) transfer step in MPM.
     *
     * @param grid_velocity The velocity field on the grid nodes.
     * @param dt Time step size.
     */
    template <typename VectorType>
    void
    grid_to_particle(const VectorType &grid_velocity, const double dt);

    /**
     * Update material point stresses based on deformation gradient.
     * This step applies the constitutive model to update stresses.
     *
     * @param dt Time step size.
     */
    void
    update_stresses(const double dt);

    /**
     * Return the number of material points.
     */
    types::particle_index
    n_global_particles() const;

    /**
     * Insert a single material point at a specified location.
     *
     * @param position The position where to insert the material point.
     * @param properties Properties of the material point (mass, volume, etc.).
     * @return Iterator to the inserted material point.
     */
    void
    insert_particle(const Point<spacedim>           &position,
                    const std::vector<double>       &properties,
                    const types::particle_index      id);

    /**
     * Update internal data structures after particle insertions.
     * This should be called after a series of insert_particle calls.
     */
    void
    update_cached_numbers();

    /**
     * Clear all material points.
     */
    void
    clear_particles();

  private:
    /**
     * Underlying particle handler for storing material points.
     */
    Particles::ParticleHandler<dim, spacedim> particle_handler;

    /**
     * DoF handler for the background grid.
     */
    DoFHandler<dim, spacedim> dof_handler;

    /**
     * Pointer to the finite element.
     */
    SmartPointer<const FiniteElement<dim, spacedim>> fe;

    /**
     * Pointer to the triangulation.
     */
    SmartPointer<const Triangulation<dim, spacedim>> triangulation;

    /**
     * Pointer to the mapping.
     */
    SmartPointer<const Mapping<dim, spacedim>> mapping;
  };

} // namespace MPM

DEAL_II_NAMESPACE_CLOSE

#endif
