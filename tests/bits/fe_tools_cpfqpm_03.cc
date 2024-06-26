// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2020 by the deal.II authors
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

#include "../tests.h"

#include "fe_tools_common.h"

// check
//   FETools::compute_projection_from_quadrature_points_matrix
// we put this into the fe_tools_common framework for simplicity, but
// in fact we ignore the second FE it passes to the check_this() function
// and we can only test as well for scalar elements, since this is all
// the function presently supports.
//
// this test simply computes the matrix and outputs some of its
// characteristics. in contrast to the _01 test, we choose a FE that has less
// DoFs than there are quadrature points -- not an uncommon case, and one that
// needs to work



template <int dim>
void
check_this(const FiniteElement<dim> &fe, const FiniteElement<dim> & /*fe2*/)
{
  // only check if both elements have
  // support points. otherwise,
  // interpolation doesn't really
  // work
  if (fe.n_components() != 1)
    return;

  // ignore this check if this FE has already
  // been treated
  static std::set<std::string> already_checked;
  if (already_checked.find(fe.get_name()) != already_checked.end())
    return;
  already_checked.insert(fe.get_name());


  // test with different quadrature formulas
  QGauss<dim> q_lhs(fe.degree + 1);
  QGauss<dim> q_rhs(fe.degree + 3);

  FullMatrix<double> X(fe.dofs_per_cell, q_rhs.size());

  FETools::compute_projection_from_quadrature_points_matrix(fe,
                                                            q_lhs,
                                                            q_rhs,
                                                            X);

  output_matrix(X);
}
