// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// like a union of the normals_* and q_point_sum_* tests: integrating
// \vec x times the normal over the surface of any body should yield
// the volume of this body times the space dimension
//
// in contrast to volume_1, do this with a Q1, Q2, and Q4 mapping
// associated with the geometry
//
// this used to fail in 3d, since there were assumptions on face
// orientations present that had to be weeded out

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



template <int dim>
void
check(const Triangulation<dim> &tria, const unsigned int order)
{
  MappingQ<dim>   mapping(order);
  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  QGauss<dim - 1> q_face(3);

  FEFaceValues<dim>    fe_face_values(mapping,
                                   fe,
                                   q_face,
                                   update_normal_vectors |
                                     update_quadrature_points |
                                     update_JxW_values);
  FESubfaceValues<dim> fe_subface_values(mapping,
                                         fe,
                                         q_face,
                                         update_normal_vectors |
                                           update_quadrature_points |
                                           update_JxW_values);

  double v1 = 0, v2 = 0;
  for (typename DoFHandler<dim>::active_cell_iterator cell =
         dof_handler.begin_active();
       cell != dof_handler.end();
       ++cell)
    {
      // first integrate over faces
      // and make sure that the
      // result of the integration is
      // close to zero
      for (const unsigned int f : GeometryInfo<dim>::face_indices())
        if (cell->at_boundary(f))
          {
            fe_face_values.reinit(cell, f);
            for (unsigned int q = 0; q < q_face.size(); ++q)
              v1 += (fe_face_values.normal_vector(q) *
                     fe_face_values.quadrature_point(q)) *
                    fe_face_values.JxW(q);
          }

      // now same for subface
      // integration
      for (const unsigned int f : GeometryInfo<dim>::face_indices())
        if (cell->at_boundary(f))
          for (unsigned int sf = 0;
               sf < GeometryInfo<dim>::max_children_per_face;
               ++sf)
            {
              fe_subface_values.reinit(cell, f, sf);
              for (unsigned int q = 0; q < q_face.size(); ++q)
                v2 += (fe_subface_values.normal_vector(q) *
                       fe_subface_values.quadrature_point(q)) *
                      fe_subface_values.JxW(q);
            }
    }

  deallog << " face integration: " << v1 / dim << std::endl;
  deallog << " subface integration: " << v2 / dim << std::endl;
  const double relative_error = std::fabs(v1 - v2) / v1;
  deallog << " relative error: " << relative_error << std::endl;
  Assert(relative_error < 5e-4, ExcInternalError());
}


int
main()
{
  initlog();

  {
    Triangulation<2> coarse_grid;
    GridGenerator::hyper_ball(coarse_grid);
    static const SphericalManifold<2> boundary;
    coarse_grid.set_manifold(0, boundary);
    check(coarse_grid, 1);
    check(coarse_grid, 2);
    check(coarse_grid, 4);
  }

  {
    Triangulation<3> coarse_grid;
    GridGenerator::hyper_ball(coarse_grid);
    static const SphericalManifold<3> boundary;
    coarse_grid.set_manifold(0, boundary);
    check(coarse_grid, 1);
    check(coarse_grid, 2);
    check(coarse_grid, 3);
  }
}
