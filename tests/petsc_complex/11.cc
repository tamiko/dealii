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



// check PETScWrappers::Vector::size()

#include "../tests.h"
#include <deal.II/lac/petsc_vector.h>
#include <fstream>
#include <iostream>


void test_real (PETScWrappers::Vector &v)
{
  deallog << "Real test" << std::endl;

  // set only certain elements of the vector
  for (unsigned int k=0; k<v.size(); k+=1+k)
    v(k) = k;

  v.compress (VectorOperation::add);

  Assert (v.size() == 100, ExcInternalError());

  deallog << "OK" << std::endl;
}

void test_complex (PETScWrappers::Vector &v)
{
  deallog << "Complex test" << std::endl;

  // set only certain elements of the vector
  for (unsigned int k=0; k<v.size(); k+=1+k)
    v(k) = std::complex<double> (k,.5*k);

  v.compress (VectorOperation::add);

  Assert (v.size() == 100, ExcInternalError());

  deallog << "OK" << std::endl;
}



int main (int argc,char **argv)
{
  std::ofstream logfile("11/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      {
        PETScWrappers::Vector v;

	v.reinit (100);
        test_real (v);

	v.reinit (100);
        test_complex (v);
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
