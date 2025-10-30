/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2025 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * ------------------------------------------------------------------------
 */

// @sect3{Include files}

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

// Include the new MPM headers
#include <deal.II/mpm/material_point.h>
#include <deal.II/mpm/mpm_handler.h>
#include <deal.II/mpm/constitutive_model.h>

#include <fstream>
#include <iostream>

// @sect3{The main program}

// This simple example demonstrates the basic usage of the MPM module.
// It creates a simple domain, initializes material points, and shows
// how to set up an MPM simulation.

using namespace dealii;

template <int dim>
class MPMExample
{
public:
  MPMExample();
  void run();

private:
  void create_grid();
  void initialize_material_points();
  void output_results(const unsigned int timestep) const;

  Triangulation<dim> triangulation;
  MappingQ<dim>      mapping;
  MPM::MPMHandler<dim> mpm_handler;
  
  const unsigned int n_timesteps = 10;
  const double       dt          = 0.01;
};


template <int dim>
MPMExample<dim>::MPMExample()
  : mapping(1)
{}


template <int dim>
void
MPMExample<dim>::create_grid()
{
  // Create a simple rectangular domain
  GridGenerator::hyper_cube(triangulation, 0.0, 1.0);
  triangulation.refine_global(3);

  std::cout << "Number of active cells: " << triangulation.n_active_cells()
            << std::endl;
}


template <int dim>
void
MPMExample<dim>::initialize_material_points()
{
  // Initialize the MPM handler with properties
  // Properties: mass, volume, velocities (dim components)
  const unsigned int n_properties = 2 + dim;
  mpm_handler.initialize(triangulation, mapping, n_properties);

  // Create a regular grid of material points in the domain
  const unsigned int n_particles_per_direction = 4;
  const double       spacing = 1.0 / (n_particles_per_direction + 1);

  types::particle_index particle_id = 0;

  if constexpr (dim == 2)
    {
      for (unsigned int i = 1; i <= n_particles_per_direction; ++i)
        for (unsigned int j = 1; j <= n_particles_per_direction; ++j)
          {
            Point<dim> position;
            position[0] = i * spacing;
            position[1] = j * spacing;

            // Properties: mass, volume, vx, vy
            std::vector<double> properties(4);
            properties[0] = 1.0;  // mass
            properties[1] = spacing * spacing; // volume (approximate)
            properties[2] = 0.0;  // velocity x
            properties[3] = 0.0;  // velocity y

            mpm_handler.insert_particle(position, properties, particle_id++);
          }
    }
  else if constexpr (dim == 3)
    {
      for (unsigned int i = 1; i <= n_particles_per_direction; ++i)
        for (unsigned int j = 1; j <= n_particles_per_direction; ++j)
          for (unsigned int k = 1; k <= n_particles_per_direction; ++k)
            {
              Point<dim> position;
              position[0] = i * spacing;
              position[1] = j * spacing;
              position[2] = k * spacing;

              // Properties: mass, volume, vx, vy, vz
              std::vector<double> properties(5);
              properties[0] = 1.0;  // mass
              properties[1] = spacing * spacing * spacing; // volume
              properties[2] = 0.0;  // velocity x
              properties[3] = 0.0;  // velocity y
              properties[4] = 0.0;  // velocity z

              mpm_handler.insert_particle(position, properties, particle_id++);
            }
    }

  // Update cached particle data
  mpm_handler.update_cached_numbers();

  std::cout << "Number of material points: "
            << mpm_handler.n_global_particles() << std::endl;
}


template <int dim>
void
MPMExample<dim>::output_results(const unsigned int timestep) const
{
  std::cout << "  Writing output for timestep " << timestep << std::endl;
  // In a full implementation, this would output particle positions
  // and properties for visualization
}


template <int dim>
void
MPMExample<dim>::run()
{
  std::cout << "Running MPM example in " << dim << "D" << std::endl;

  create_grid();
  initialize_material_points();

  // Time stepping loop
  for (unsigned int timestep = 0; timestep < n_timesteps; ++timestep)
    {
      std::cout << "Timestep " << timestep << " at time "
                << timestep * dt << std::endl;

      // In a full MPM simulation, each timestep would include:
      // 1. Particle to grid transfer (P2G)
      // 2. Solve equations of motion on grid
      // 3. Update grid velocities
      // 4. Grid to particle transfer (G2P)
      // 5. Update particle stresses
      
      output_results(timestep);
    }

  std::cout << "MPM simulation completed." << std::endl;
}


int main()
{
  try
    {
      MPMExample<2> mpm_example_2d;
      mpm_example_2d.run();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
