// Copyright (C) 2004-2014  Carlo de Falco
// Copyright (C) 2014  Davide Cagnoni
//
// This file is part of:
//     secs3d - A 3-D Drift--Diffusion Semiconductor Device Simulator
//
//  secs3d is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  secs3d is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with secs3d; If not, see <http://www.gnu.org/licenses/>.
//
//  author: Carlo de Falco     <cdf _AT_ users.sourceforge.net>

/*! \file operators.h
  \brief Functions to build the discrete version of differential operators.
*/

#ifndef HAVE_OPERATORS_H
#define HAVE_OPERATORS_H 1

#include "mesh.h"
#include "bim_sparse.h"
#include <cmath>
#include <mpi.h>

//namespace bim
//{

/// Virtual class defining the interface
/// of a FEM coefficient functor.
class
coefficient_functor
{
public:
  virtual const double& operator () (const int) = 0;
};

/// Specialization to use a vector as a
/// coefficient functor.
class
coefficient_vector : public coefficient_functor
{
private:
  const std::vector<double> *v;
public:
  /// Constructor
  coefficient_vector (const std::vector<double> *_v) : v (_v) {};
  ///
  const double& operator () (const int ii) { return (*v)[ii]; };
};

/// Enumeration to specificate type of norm
enum norm_type {Inf,L2,H1};

/// Allocate the structure for a FEM matrix over the mesh msh.
void
bim3a_structure (const mesh& msh,
                 sparse_matrix& A);

/// Assemble the rhs of a FEM problem.
void
bim3a_rhs (mesh& msh,
           const std::vector<double>& ecoeff,
           const std::vector<double>& ncoeff,
           std::vector<double>& b);

/// Assemble the mass matrix of a FEM problem.
void
bim3a_reaction (mesh& msh,
                const std::vector<double>& ecoeff,
                const std::vector<double>& ncoeff,
                sparse_matrix& A);

/// Assemble the stiffness matrix of a FEM diffusion problem.
void
bim3a_laplacian (mesh& msh,
                 const std::vector<double>& epsilon,
                 sparse_matrix& A);

/// Assemble the stiffness matrix of a FEM diffusion problem
/// with anisotropic diffusion.
void
bim3a_laplacian_anisotropic (mesh& msh,
                             const std::vector<double>& epsilon,
                             sparse_matrix& A);

/// Assemble the stiffness matrix of a FEM advection-diffusion problem.
void
bim3a_advection_diffusion (mesh& msh,
                           const std::vector<double>& epsilon,
                           const std::vector<double>& phi,
                           sparse_matrix& A);

/// Assemble the stiffness matrix of a FEM advection-diffusion problem
/// with anisotropic diffusion.
void
bim3a_advection_diffusion_anisotropic (mesh& msh,
                                       const std::vector<double>& epsilon,
                                       const std::vector<double>& phi,
                                       sparse_matrix& A);

/// Assemple the stiffness matrix of a FEM advection problem
/// with upwind stabilisation.
void
bim3a_advection_upwind (mesh &msh,
                        const std::vector<double>& v,
                        sparse_matrix& UP);

/// Assemble the stiffness matrix of a FEM diffusion problem
/// (by the Orthogonal Subdomain Collocation method).
void
bim3a_osc_laplacian (mesh& msh,
                     const std::vector<double>& epsilon,
                     sparse_matrix& A);

/// Assemble the stiffness matrix of a FEM diffusion problem
/// with anisotropic diffusion
/// (by the Orthogonal Subdomain Collocation method).
void
bim3a_osc_laplacian_anisotropic (mesh& msh,
                                 const std::vector<double>& epsilon,
                                 sparse_matrix& A);

/// Assemble the stiffness matrix of a FEM advection-diffusion
/// problem (by the Orthogonal Subdomain Collocation method).
void
bim3a_osc_advection_diffusion (mesh& msh,
                               const std::vector<double>& epsilon,
                               const std::vector<double>& phi,
                               sparse_matrix& A);

/// Assemble the stiffness matrix of a FEM advection-diffusion
/// problem with anisotropic diffusion
/// (by the Orthogonal Subdomain Collocation method).
void
bim3a_osc_advection_diffusion_anisotropic (mesh& msh,
                                           const std::vector<double>& epsilon,
                                           const std::vector<double>& phi,
                                           sparse_matrix& A);

/// Robustly compute B(x) = x / (exp(x) - 1). Stores B(x) and B(-x)
void
bimu_bernoulli (double x,
                double &bp,
                double &bm);

/// Robustly compute B'(x), B(x) = x / (exp(x) - 1).
/// Stores B'(x) and B'(-x)
void
bimu_bernoulli_derivative (double x,
                           double &bpp,
                           double &bmp);

/// Modify a precomputed local diffusion matrix, considering
/// effects due to advection field generated by a potential,
/// through Scharfetter-Gummel method
void
bim3a_local_advection (const double v[4],
                       const double a,
                       const double epsilon,
                       double Lloc[16]);

/// Transform a precomputed local diffusion matrix into the jacobian
/// of the S-G transport term with respect to the potential of the
/// transport field
void
bim3a_local_advection_jacobian (const double v[4],
                                const double n[4],
                                const double a,
                                const double epsilon,
                                double Lloc[16]);

/// Compute the elemental contribution to the global
/// stiffness matrix (FEM), and add it to the exit buffer
void
bim3a_local_laplacian (const double shg[12],
                       const double volume,
                       const double epsilon,
                       double Lloc[16]);

/// Compute the elemental contribution to the global
/// stiffness matrix (FEM) with anisotropic diffusion,
/// and add it to the exit buffer
void
bim3a_local_laplacian_anisotropic (const double shg[12],
                                   const double volume,
                                   const double epsilon[9],
                                   double Lloc[16]);

/// Compute the elemental contribution to the global
/// stiffness matrix (FEM) with diagonal anisotropic diffusion,
/// and add it to the exit buffer
void
bim3a_local_laplacian_anisotropic_diag (const double shg[12],
                                        const double volume,
                                        const double epsilon[3],
                                        double Lloc[16]);

/// Compute the elemental contribution to the global
/// stiffness matrix (OSC method), and add it to the exit buffer
void
bim3a_osc_local_laplacian (const double shg[12],
                           const double p[12],
                           const double volume,
                           const double epsilon,
                           double Lloc[16]);

/// Compute the elemental contribution to the global
/// stiffness matrix (OSC method) with diagonal anisotropic diffusion,
/// and add it to the exit buffer
void
bim3a_osc_local_laplacian_anisotropic (const double shg[12],
                                       const double p[12],
                                       const double volume,
                                       const double epsilon[3],
                                       double Lloc[16]);

/// Compute the elemental contribution to the global mass matrix
/// (with mass lumping), and add it to the exit buffer
void
bim3a_local_reaction (const double shp[16],
                      const double wjacdet[4],
                      const double coeffe,
                      const double coeffn[4],
                      double Lloc[16]);

/// Compute the elemental contribution to the rhs of a FEM problem
/// and add it to the exit buffer.
void
bim3a_local_rhs (const double shp[16],
                 const double wjacdet[4],
                 const double coeffe,
                 const double coeffn[4],
                 double bLoc[4]);

/// Set Dirichlet boundary condition.
void
bim3a_dirichlet_bc (sparse_matrix& M,
                    std::vector<double>& b,
                    const std::vector<int>& bnodes,
                    const std::vector<double>& vnodes);

/// Compute the gradient of a piecewise linear function.
void
bim3a_pde_gradient (mesh& msh,
                    const std::vector<double>& u,
                    std::vector<double>& g);

/// Compute the (Inf,L2,H1) norm of a piecewise linear function.
void
bim3a_norm (mesh& msh,
            const std::vector<double>& v,
            double& norm,
            norm_type type);

//}
#endif
