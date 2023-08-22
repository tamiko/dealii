// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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

// Extract cell bounding boxes rtree from the cache, and try to use it

#include <deal.II/base/patterns.h>

#include <deal.II/boost_adaptors/point.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools_cache.h>
#include <deal.II/grid/tria.h>

#include <algorithm>

#include "../tests.h"

using Patterns::Tools::to_string;

namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;

template <int dim, int spacedim>
void
test(const unsigned int ref = 2, const unsigned int n_points = 10)
{
  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(ref);
  GridTools::Cache<dim, spacedim> cache(tria);

  const auto &b_tree = cache.get_cell_bounding_boxes_rtree();

  std::vector<Point<spacedim>> points(n_points);
  std::generate(points.begin(), points.end(), []() { return random_point<spacedim>(); });

  deallog << "Testing dim = " << dim << ", spacedim = " << spacedim << std::endl;

  for (const auto &p : points)
    {
      std::vector<std::pair<BoundingBox<spacedim>, typename Triangulation<dim, spacedim>::active_cell_iterator>> res;
      b_tree.query(bgi::nearest(p, 1), std::back_inserter(res));
      deallog << "Nearest cell to " << p << ":  " << res[0].second << std::endl;
    }
}

int
main(int argc, char **argv)
{
  initlog();

  test<1, 1>();
  test<1, 2>();
  test<1, 3>();
  test<2, 2>();
  test<2, 3>();
  test<3, 3>();
}
