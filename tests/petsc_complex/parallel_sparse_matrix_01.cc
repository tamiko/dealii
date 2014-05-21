// ---------------------------------------------------------------------
// $Id: testsuite.html 32446 2014-02-10 14:31:22Z bangerth $
//
// Copyright (C) 2013 by the deal.II Authors
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

// test that matrix_scalar_product of a symmetric matrix 
// applied to the same vector result in a real number


#include "../tests.h"

#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>            
#include <deal.II/base/conditional_ostream.h>  
#include <deal.II/base/index_set.h>            
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>          
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>

#include <deal.II/distributed/tria.h>            
#include <deal.II/distributed/grid_refinement.h> 

using namespace dealii;

static const unsigned int dim = 2;

void test ()
{

  MPI_Comm mpi_communicator (MPI_COMM_WORLD);
   const unsigned int n_mpi_processes (Utilities::MPI::n_mpi_processes(mpi_communicator));
   const unsigned int this_mpi_process (Utilities::MPI::this_mpi_process(mpi_communicator));
   
   ConditionalOStream pcout (std::cout,
           (Utilities::MPI::this_mpi_process(mpi_communicator)
            == 0));
  
  parallel::distributed::Triangulation<dim> tria(mpi_communicator,
                   typename Triangulation<dim>::MeshSmoothing
                   (Triangulation<dim>::smoothing_on_refinement |
                    Triangulation<dim>::smoothing_on_coarsening));

  if (false)
   {
    Triangulation<dim>   triangulation1;
    Triangulation<dim>   triangulation2;
    Triangulation<dim>   triangulation12;
    GridGenerator::hyper_cube (triangulation1, -1,0); //create a square [-1,0]^d domain
    GridGenerator::hyper_cube (triangulation2, -1,0); //create a square [-1,0]^d domain
    Point<dim> shift_vector;
    shift_vector[0] = 1.0;
    GridTools::shift(shift_vector,triangulation2);
    GridGenerator::merge_triangulations (triangulation1, triangulation2, tria);
   }  
  else
  {
     GridGenerator::hyper_cube (tria, -1,0);
	 //tria.refine_global (1);
  } 
  
  const unsigned int poly_degree = 1;
  FE_Q<dim> fe(poly_degree);
      
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);
  
  IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs ();
  IndexSet locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs (dof_handler,locally_relevant_dofs);
  
  PETScWrappers::MPI::Vector vector, conjugate;
  PETScWrappers::MPI::SparseMatrix mass_matrix;
  

  vector.reinit(locally_owned_dofs, mpi_communicator);
  const types::global_dof_index n_local_dofs = locally_owned_dofs.n_elements();
  mass_matrix.reinit  (mpi_communicator,
			   dof_handler.n_dofs(),
			   dof_handler.n_dofs(),
			   n_local_dofs,
			   n_local_dofs,
			   dof_handler.max_couplings_between_dofs());
			   
  ConstraintMatrix constraints;			   
  constraints.reinit(locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints  (dof_handler, constraints);
  constraints.close ();
			   
  //assemble mass matrix:
  mass_matrix = 0.0;
	  {
		  //a separate quadrature formula
		  //enough for mass and kinetic matrices assembly
		  QGauss<dim> quadrature_formula(poly_degree+1);
		  FEValues<dim> fe_values (fe, quadrature_formula,
				                   update_values | update_gradients | update_JxW_values);

		  const unsigned int dofs_per_cell = fe.dofs_per_cell;
		  const unsigned int n_q_points    = quadrature_formula.size();

		  FullMatrix<PetscScalar> cell_mass_matrix (dofs_per_cell, dofs_per_cell);

		  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

		  typename DoFHandler<dim>::active_cell_iterator
			  cell = dof_handler.begin_active (),
			  endc = dof_handler.end ();
		  for (; cell!=endc; ++cell)
			if (cell->is_locally_owned())
			{
			   fe_values.reinit (cell);
			   cell_mass_matrix    = PetscScalar(0);

			   for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
				 for (unsigned int i=0; i<dofs_per_cell; ++i)
				   for (unsigned int j=0; j<dofs_per_cell; ++j)
					 {
					   cell_mass_matrix (i, j)
					   += (fe_values.shape_value (i, q_point) *
						   fe_values.shape_value (j, q_point)
						  ) * fe_values.JxW (q_point);
					 }

			   cell->get_dof_indices (local_dof_indices);

			   constraints
			   .distribute_local_to_global (cell_mass_matrix,
											local_dof_indices,
											mass_matrix);
			 }
			 mass_matrix.compress (VectorOperation::add);
	  }			   
  
  vector = std::complex<double>(1.0,-1.0);
  
  //deallog<<"original vector: "<<std::endl;
  //vector.print(std::cout);
  
  deallog<<"symmetric: "<<mass_matrix.is_symmetric()<<std::endl;
  deallog<<"hermitian: "<<mass_matrix.is_hermitian()<<std::endl;
  
  //mass_matrix.print(std::cout);
  
  PETScWrappers::MPI::Vector tmp(vector);
  mass_matrix.vmult (tmp, vector);
  
  //pcout<<"vmult: "<<std::endl;
  //tmp.print(std::cout);
  //correct output ^^^ is:
  // 0.25000 - 0.25000i
  // 0.25000 - 0.25000i
  // 0.25000 - 0.25000i
  // 0.25000 - 0.25000i
  
  conjugate = vector;
  conjugate.conjugate();
  //pcout<<"conjugate vector: "<<std::endl;
  //conjugate.print(std::cout);
  
  
  deallog<<"conj*vmult:"<< conjugate*tmp<<std::endl;
  
  deallog<<"vector*vmult:"<< vector*tmp<<std::endl;
  
  deallog<<"matrix_scalar_product: ";
  const std::complex<double> norm =
	  mass_matrix.matrix_scalar_product(vector, vector);
  deallog<<norm<<std::endl;
  //correct output is:
  // (2,0)
}//*/


int main (int argc, char *argv[])
{
   std::ofstream logfile("output");
   deallog.attach(logfile);
   deallog.depth_console(0);
   deallog.threshold_double(1.e-10);
  
   Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv);
   
   test ();
   
}