# Material Point Method (MPM) Module for deal.II

This module implements the Material Point Method (MPM) for deal.II, a computational technique for simulating solid mechanics problems involving large deformations, contact, and material failure.

## Overview

The Material Point Method is a particle-based method that combines the advantages of Lagrangian and Eulerian approaches. Material points carry all the material state (mass, volume, velocity, stress, strain), while a background Eulerian grid is used to solve the equations of motion.

## Key Components

### Material Point (`material_point.h`)
Represents a single material point with properties:
- Mass
- Volume  
- Velocity
- Stress tensor
- Strain tensor

Extends the base `Particles::Particle` class with material-specific data.

### MPM Handler (`mpm_handler.h`)
Manages collections of material points and provides core MPM operations:
- `particle_to_grid()` - Transfer mass and momentum from material points to grid (P2G)
- `compute_grid_forces()` - Compute internal forces from material point stresses
- `grid_to_particle()` - Update material points from grid solution (G2P)
- `update_stresses()` - Update material point stresses using constitutive models

### Constitutive Model (`constitutive_model.h`)
Material constitutive models that relate stress to strain:
- `ConstitutiveModel` - Abstract base class for constitutive laws
- `LinearElasticModel` - Linear elastic material model

## Example Usage

See `examples/step-mpm/step-mpm.cc` for a basic example that demonstrates:
1. Creating a computational domain
2. Initializing material points
3. Setting up an MPM simulation framework

## Implementation Status

This is an initial implementation providing the basic framework for MPM simulations. The current implementation includes:

✅ Material point data structure
✅ MPM handler for managing material points
✅ Linear elastic constitutive model
✅ Framework for P2G and G2P transfers (placeholders)
✅ Basic example

Future enhancements could include:
- Complete implementation of P2G and G2P transfer operations
- Additional constitutive models (plasticity, hyperelasticity)
- Contact algorithms
- Damage and fracture models
- Parallel implementation using MPI

## References

1. Sulsky, D., Chen, Z., & Schreyer, H. L. (1994). A particle method for history-dependent materials. Computer methods in applied mechanics and engineering, 118(1-2), 179-196.

2. Bardenhagen, S. G., & Kober, E. M. (2004). The generalized interpolation material point method. Computer Modeling in Engineering and Sciences, 5(6), 477-496.

3. Stomakhin, A., Schroeder, C., Chai, L., Teran, J., & Sifakis, E. (2013). A material point method for snow simulation. ACM Transactions on Graphics (TOG), 32(4), 1-10.
