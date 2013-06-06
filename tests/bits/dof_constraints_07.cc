//----------------------------  dof_constraints_07.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  dof_constraints_07.cc  ---------------------------


// simply check what happens when calling DoFConstraints::set_zero on
// block vectors. This test was written when I changed a few things in the algorithm

#include "../tests.h"
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/fe/fe_q.h>
#include <fstream>
#include <iomanip>


template <int dim>
void test ()
{
  deallog << dim << "D" << std::endl;
  
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube (triangulation);

                                   // refine once, then refine first cell to
                                   // create hanging nodes
  triangulation.refine_global (1);
  triangulation.begin_active()->set_refine_flag ();
  triangulation.execute_coarsening_and_refinement ();
  deallog << "Number of cells: " << triangulation.n_active_cells() << std::endl;
  
                                   // set up a DoFHandler and compute hanging
                                   // node constraints for a Q2 element
  FE_Q<dim> fe(2);
  DoFHandler<dim> dof_handler (triangulation);
  dof_handler.distribute_dofs (fe);
  deallog << "Number of dofs: " << dof_handler.n_dofs() << std::endl;

  ConstraintMatrix constraints;
  DoFTools::make_hanging_node_constraints (dof_handler, constraints);
  constraints.close ();
  deallog << "Number of constraints: " << constraints.n_constraints() << std::endl;

  std::vector<types::global_dof_index> block_sizes(2);
  block_sizes[0] = dof_handler.n_dofs()/3;
  block_sizes[1] = dof_handler.n_dofs() - block_sizes[0];
  BlockVector<double> b(block_sizes);
  for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
    b(i) = (1.+1.*i*i)/3;

                                   // now condense away constraints
  constraints.set_zero (b);

                                   // and output what we have
  for (BlockVector<double>::const_iterator i=b.begin(); i!=b.end(); ++i)
    deallog << *i << std::endl;

                                   // now also make sure that the elements in
                                   // constrained rows are zero, and that all
                                   // the other elements are unchanged
  for (unsigned int i=0; i<b.size(); ++i)
    if (constraints.is_constrained(i))
      Assert (b(i) == 0, ExcInternalError())
    else
      Assert (std::fabs(b(i) - (1.+1.*i*i)/3) < 1e-14*std::fabs(b(i)),
              ExcInternalError());
}



int main ()
{
  std::ofstream logfile("dof_constraints_07/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  try
    {
      test<1> ();
      test<2> ();
      test<3> ();
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
