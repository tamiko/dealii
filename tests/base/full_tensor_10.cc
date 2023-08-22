// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2018 by the deal.II authors
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


// test the determinant code for n>3

#include <deal.II/base/tensor.h>

#include "../tests.h"


template <int dim>
void
test()
{
  Tensor<2, dim> t;

  // choose the same symmetric tensor
  // as in symmetric_tensor_10
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      t[i][j] = random_value<double>();

  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      deallog << "A[" << i + 1 << ',' << j + 1 << "] := " << t[i][j] << ';' << std::endl;

  deallog << determinant(t) << std::endl;
}



int
main()
{
  initlog();
  deallog << std::setprecision(3);

  test<4>();
  test<5>();
  test<6>();
  test<7>();
  test<8>();

  deallog << "OK" << std::endl;
}
