// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2018 by the deal.II authors
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



// check querying the number of nonzero elements in
// SparseMatrix when we don't store the diagonal elements explicitly

#include <deal.II/lac/sparse_matrix.h>

#include "../tests.h"


void
test()
{
  SparsityPattern sp(5, 5, 3);
  for (unsigned int i = 0; i < 5; ++i)
    for (unsigned int j = 0; j < 5; ++j)
      if ((i + 2 * j + 1) % 3 == 0)
        sp.add(i, j);
  sp.compress();

  SparseMatrix<double> m(sp);

  // first set a few entries. count how many
  // entries we have. note that for square
  // matrices we also always store the
  // diagonal element
  unsigned int counter = 0;
  for (unsigned int i = 0; i < m.m(); ++i)
    for (unsigned int j = 0; j < m.n(); ++j)
      if (((i + 2 * j + 1) % 3 == 0) || (i == j))
        {
          m.set(i, j, i * j * .5 + .5);
          ++counter;
        }

  deallog << m.n_nonzero_elements() << std::endl;
  Assert(m.n_nonzero_elements() == counter, ExcInternalError());

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  try
    {
      test();
    }
  catch (const std::exception &exc)
    {
      deallog << std::endl << std::endl << "----------------------------------------------------" << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------" << std::endl;

      return 1;
    }
  catch (...)
    {
      deallog << std::endl << std::endl << "----------------------------------------------------" << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------" << std::endl;
      return 1;
    };
}
