//----------------------------  trilinos_67.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2004, 2005, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  trilinos_67.cc  ---------------------------


// check TrilinosWrappers::MatrixBase::clear_rows ()

#include "../tests.h" 
#include <deal.II/base/utilities.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <fstream>
#include <iostream>
#include <vector>


void test (TrilinosWrappers::SparseMatrix &m)
{
  Assert (m.m() != 0, ExcInternalError());
  Assert (m.n() != 0, ExcInternalError());

                                   // build a tri-diagonal pattern
  double norm_sqr = 0;
  unsigned int nnz = 0;
  const unsigned int N = m.m();
  for (unsigned int i=0; i<N; ++i)
    {
      if (i>=5)
        {
          const double s = rand();
          m.set (i,i-5, s);
          norm_sqr += s*s;
          ++nnz;
        }
      
      if (i<N-5)
        {
          const double s = rand();
          m.set (i,i+5, s);
          norm_sqr += s*s;
          ++nnz;
        }
      
      const double s = rand();
      m.set (i,i,s);
      norm_sqr += s*s;
      ++nnz;
    }
  m.compress (VectorOperation::insert);
  
  deallog << m.frobenius_norm() << ' ' << std::sqrt (norm_sqr)
          << std::endl;
  deallog << m.n_nonzero_elements() << ' ' << nnz << std::endl;

  Assert (std::fabs (m.frobenius_norm() - std::sqrt(norm_sqr))
          < std::fabs (std::sqrt (norm_sqr)),
          ExcInternalError());
  Assert (m.n_nonzero_elements()-nnz == 0, ExcInternalError());

                                   // now remove the entries of rows N/2 and N/3
  for (unsigned int i=0; i<N; ++i)
    {
      const double s = m.el(N/2,i);
      norm_sqr -= s*s;
    }
  for (unsigned int i=0; i<N; ++i)
    {
      const double s = m.el(N/3,i);
      norm_sqr -= s*s;
    }
  const types::global_dof_index rows[2] = { N/3, N/2 };
  m.clear_rows (std::vector<types::global_dof_index>(&rows[0], &rows[2]));
  
  deallog << m.frobenius_norm() << ' ' << std::sqrt (norm_sqr)
          << std::endl;
  deallog << m.n_nonzero_elements() << ' ' << nnz << std::endl;

  Assert (std::fabs (m.frobenius_norm() - std::sqrt(norm_sqr))
          < std::fabs (std::sqrt (norm_sqr)),
          ExcInternalError());

                                   // make sure that zeroing out rows does at
                                   // least not add new nonzero entries (it
                                   // may remove some, though)
  Assert (m.n_nonzero_elements() <= nnz, ExcInternalError());
  
  deallog << "OK" << std::endl;
}



int main (int argc,char **argv) 
{
  std::ofstream logfile("67/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10); 

  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv);


  try
    {
      {
        TrilinosWrappers::SparseMatrix v (14U,14U,3U);
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
