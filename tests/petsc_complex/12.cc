// ---------------------------------------------------------------------
// $Id: 12.cc $
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



// check PETScWrappers::Vector::operator() in set-mode

#include "../tests.h"
#include <deal.II/lac/petsc_vector.h>
#include <fstream>
#include <iostream>
#include <vector>



void test (PETScWrappers::Vector &v)
{
  deallog << "Complex test" << std::endl;

  // set only certain elements of the vector. have a bit pattern of
  // where we actually wrote elements to
  std::vector<bool> pattern (v.size(), false);
  for (unsigned int k=0; k<v.size(); k+=1+k)
    {
      v(k)       = std::complex<double> (0.,k);
      pattern[k] = true;
    }

  v.compress (VectorOperation::add);

  // check that they are ok, and this time all of them
  for (unsigned int k=0; k<v.size(); ++k)
    {
      const PetscScalar el = v(k);
      Assert ( ( (pattern[k] == true) && (PetscRealPart(el) == 0.) && (PetscImaginaryPart(el) == k) )
	       ||
	       ( (pattern[k] == false) && (PetscRealPart(el) == 0.) && (PetscImaginaryPart(el) == 0.)),
	       ExcInternalError());
    }

  deallog << "OK" << std::endl;
}


int main (int argc,char **argv)
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      {
        PETScWrappers::Vector v(100);
        test (v);
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
