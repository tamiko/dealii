// ---------------------------------------------------------------------
// $Id: 
//
// Copyright (C) 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// check querying the number of nonzero elements in
// PETScWrappers::SparseMatrix

#include "../tests.h"
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <fstream>
#include <iostream>


void test_real (PETScWrappers::SparseMatrix &m)
{
  deallog << "Real test" << std::endl;

  // first set a few entries. count how many entries we have
  unsigned int counter = 0;
  for (unsigned int i=0; i<m.m(); ++i)
    for (unsigned int j=0; j<m.m(); ++j)
      if ((i+2*j+1) % 3 == 0)
        {
          m.set (i,j, i*j*.5+.5);
          ++counter;
        }

  m.compress (VectorOperation::add);

  deallog << m.n_nonzero_elements() << std::endl;
  Assert (m.n_nonzero_elements() == counter,
          ExcInternalError());

  deallog << "OK" << std::endl;
}


void test_complex (PETScWrappers::SparseMatrix &m)
{
  deallog << "Complex test" << std::endl;

  // first set a few entries. count how many entries we have
  unsigned int counter = 0;
  for (unsigned int i=0; i<m.m(); ++i)
    for (unsigned int j=0; j<m.m(); ++j)
      if ((i+2*j+1) % 3 == 0)
        {
          m.set (i,j, std::complex<double> (0.,i*j*.5+.5));
          ++counter;
        }

  m.compress (VectorOperation::add);

  deallog << m.n_nonzero_elements() << std::endl;
  Assert (m.n_nonzero_elements() == counter,
          ExcInternalError());

  deallog << "OK" << std::endl;
}



int main (int argc,char **argv)
{
  std::ofstream logfile("05/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      {
        PETScWrappers::SparseMatrix m (5,5,3);
        test_real (m);

        PETScWrappers::SparseMatrix n (5,5,3);
        test_complex (n);

	// this last bit is unnecessarry fun...
	deallog << "Compare real-complex test" << std::endl;
	Assert (m.n_nonzero_elements() == n.n_nonzero_elements(),
		ExcInternalError());
	deallog << "OK" << std::endl;
      }

    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
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
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    };
}
