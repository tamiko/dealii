//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


// create a parallel DoFHandler on a single CPU

#include "../tests.h"
#include "coarse_grid_common.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_q.h>

#include <fstream>


template<int dim>
void test(std::ostream& /*out*/)
{
  deallog << "hyper_cube" << std::endl;

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);
      
  GridGenerator::hyper_cube(tr);
  tr.refine_global (1);
  DoFHandler<dim> dofh(tr);

  static const FE_Q<dim> fe(2);
  dofh.distribute_dofs (fe);

   typename
     DoFHandler<dim>::active_cell_iterator cell
     = dofh.begin_active();

   const unsigned int dofs_per_cell = dofh.get_fe().dofs_per_cell;
   std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);


      
   for (; cell != dofh.end(); ++cell)
     {
       cell->get_dof_indices (local_dof_indices);

       for (unsigned int i=0;i<dofs_per_cell;++i)
 	deallog << local_dof_indices[i] << " ";

       deallog << std::endl;
     }
  
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  
  std::ofstream logfile("2d_dofhandler_01/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  deallog.push("2d");
  test<2>(logfile);
  deallog.pop();


}
