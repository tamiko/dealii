// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2018 by the deal.II authors
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


// check GridTools::volume for codim-one


#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

template <int dim>
void
test1()
{
  // test 1: hypercube
  if (true)
    {
      Triangulation<dim, dim + 1> tria;
      GridGenerator::hyper_cube(tria);

      for (unsigned int i = 0; i < 2; ++i)
        {
          tria.refine_global(2);
          deallog << dim << "d, "
                  << "hypercube volume, " << i * 2 << " refinements: " << GridTools::volume(tria) << std::endl;
        };
    };
}


int
main()
{
  initlog();
  deallog << std::setprecision(4);

  test1<1>();
  test1<2>();

  return 0;
}
