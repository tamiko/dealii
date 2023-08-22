// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



// this function tests the correctness of the implementation of matrix free
// matrix-vector products by comparing with the result of deal.II sparse
// matrix. The mesh uses a hyperball mesh with hanging nodes generated by
// randomly refining some cells. Similar to matrix_vector_06, but uses FE_DGQ
// which do not have any connections because DoFs do not overlap and we do not
// evaluate flux terms (this stresses the partitioning in a different way
// where non-connected DoFs are to be handled)

#include <deal.II/base/function.h>

#include <deal.II/fe/fe_dgq.h>

#include "../tests.h"

#include "matrix_vector_common.h"


template <int dim, int fe_degree>
void
test()
{
  const SphericalManifold<dim> manifold;
  Triangulation<dim>           tria;
  GridGenerator::hyper_ball(tria);
  typename Triangulation<dim>::active_cell_iterator cell = tria.begin_active(), endc = tria.end();
  for (; cell != endc; ++cell)
    for (const unsigned int f : GeometryInfo<dim>::face_indices())
      if (cell->at_boundary(f))
        cell->face(f)->set_all_manifold_ids(0);
  tria.set_manifold(0, manifold);
  cell = tria.begin_active();
  for (; cell != endc; ++cell)
    if (cell->center().norm() < 1e-8)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  cell = tria.begin_active();
  for (; cell != endc; ++cell)
    if (cell->center().norm() < 0.2)
      cell->set_refine_flag();
  tria.execute_coarsening_and_refinement();
#ifndef DEBUG
  if (dim < 3 || fe_degree < 2)
    tria.refine_global(1);
  tria.begin(tria.n_levels() - 1)->set_refine_flag();
  tria.last()->set_refine_flag();
  tria.execute_coarsening_and_refinement();
  tria.refine_global(4 - dim);
#endif
  cell = tria.begin_active();
  for (unsigned int i = 0; i < 10 - 3 * dim; ++i)
    {
      cell                 = tria.begin_active();
      unsigned int counter = 0;
      for (; cell != endc; ++cell, ++counter)
        if (counter % (7 - i) == 0)
          cell->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }

  FE_DGQ<dim>     fe(fe_degree);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);
  AffineConstraints<double> constraints;
  // there should not be any hanging nodes or
  // boundary conditions for FE_DGQ as there are
  // only interior DoFs on the elements, but try
  // anyway
  DoFTools::make_hanging_node_constraints(dof, constraints);
  VectorTools::interpolate_boundary_values(dof, 0, Functions::ZeroFunction<dim>(), constraints);
  constraints.close();

  do_test<dim, fe_degree, double, fe_degree + 1>(dof, constraints);

  // test with coloring only as well
  do_test<dim, fe_degree, double, fe_degree + 1>(dof, constraints, 2);
}
