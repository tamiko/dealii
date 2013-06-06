//----------------------------  sparse_matrix_entries_04.cc  -----------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  sparse_matrix_add_entries_04.cc  -------------


// check adding elements into a matrix using
// SparseMatrix::add(row, n_cols, col_indices, values, elide_zero_values,
//                   col_indices_are_sorted)
// rectangular case, sorted

#include "../tests.h"
#include <deal.II/lac/sparse_matrix.h>
#include <fstream>


void test ()
{
				// set up sparse matrix
  SparsityPattern sp (5,7,3);
  for (unsigned int i=0; i<sp.n_rows(); ++i)
    for (unsigned int j=0; j<sp.n_cols(); ++j)
      if ((i+2*j+1) % 3 == 0)
        sp.add (i,j);
  sp.compress ();

  SparseMatrix<double> m(sp);

				// prepare structure with indices and values
  std::vector<types::global_dof_index> indices (m.n());
  for (unsigned int j=0; j<m.n(); ++j)
    indices[j] = j;
  std::vector<double> values (m.n());

                                   // try to add entries from the list. Zeros
                                   // should be filtered out. list is sorted
  for (unsigned int i=0; i<m.m(); ++i)
    {
      for (unsigned int j=0; j<m.n(); ++j)
	if ((i+2*j+1) % 3 == 0)
	  values[j] = i*j*.5+.5;
	else
	  values[j] = 0;
      m.add(i,m.n(),&indices[0], &values[0], false, true);
    }

                                   // then make sure we retrieve the same ones
  for (unsigned int i=0; i<m.m(); ++i)
    for (unsigned int j=0; j<m.n(); ++j)
      if ((i+2*j+1) % 3 == 0)
	{
	  Assert (m(i,j) == i*j*.5+.5, ExcInternalError());
	}
      else
	{
	  Assert (m.el(i,j) == 0, ExcInternalError());
	}

  deallog << "OK" << std::endl;
}



int main ()
{
  std::ofstream logfile("sparse_matrix_add_entries_04/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      test ();
    }
  catch (std::exception &exc)
    {
      deallog << std::endl << std::endl
	      << "----------------------------------------------------"
	      << std::endl;
      deallog << "Exception on processing: " << std::endl
	      << exc.what() << std::endl
	      << "Aborting!" << std::endl
	      << "----------------------------------------------------"
	      << std::endl;

      return 1;
    }
  catch (...)
    {
      deallog << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      deallog << "Unknown exception!" << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    };
}
