//-----------------  sparse_matrix_copy_from_01.cc  -------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2011, 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------  sparse_matrix_copy_from_01.cc  -------------------------


// test SparseMatrix::copy_from from a TrilinosWrappers::SparseMatrix

#include "../tests.h"
#include <deal.II/base/utilities.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <fstream>
#include <iostream>


int main (int argc,char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv);

  std::ofstream logfile("sparse_matrix_copy_from_01/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  SparsityPattern sparsity (5,5,5);
  sparsity.add (1,2);
  sparsity.add (2,3);
  sparsity.add (3,2);
  sparsity.add (3,3);
  sparsity.add (3,4);
  sparsity.add (4,3);
  sparsity.compress();
  SparseMatrix<double> matrix(sparsity);
  {
    double value = 1;
    for (SparseMatrix<double>::iterator p=matrix.begin();
	 p != matrix.end(); ++p, ++value)
      p->value() = value;
  }
  deallog << "Original:" << std::endl;
  matrix.print_formatted (deallog.get_file_stream());

  // now copy everything into a Trilinos matrix
  Epetra_Map map(TrilinosWrappers::types::int_type(5),5,0,Utilities::Trilinos::comm_world());
  TrilinosWrappers::SparseMatrix tmatrix;
  tmatrix.reinit (map, map, matrix);

  // now copy things back into a SparseMatrix
  SparseMatrix<double> copy (sparsity);
  copy.copy_from (tmatrix);

  // print some output
  deallog << "Copy with all values:" << std::endl;
  matrix.print (deallog.get_file_stream());

  // also compare for equality with the original
  for (SparsityPattern::const_iterator
	 p = sparsity.begin(); p != sparsity.end(); ++p)
    Assert (copy(p->row(), p->column()) == matrix(p->row(), p->column()),
	    ExcInternalError());
}
