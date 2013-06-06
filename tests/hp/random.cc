//----------------------------  random.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  random.cc  ---------------------------


// see what happens if we distribute different finite elements
// randomly across the domain


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_q.h>

#include <fstream>



template <int dim>
void test ()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global (2);
  tria.begin_active()->set_refine_flag ();
  tria.execute_coarsening_and_refinement ();
  tria.refine_global (2);  

  hp::FECollection<dim> fe_collection;
  fe_collection.push_back(FE_Q<dim> (1));
  fe_collection.push_back(FE_Q<dim> (2));
  fe_collection.push_back(FE_Q<dim> (3));
  fe_collection.push_back(FE_Q<dim> (4));

  hp::DoFHandler<dim> dof_handler(tria);

  for (typename hp::DoFHandler<dim>::active_cell_iterator
	 cell = dof_handler.begin_active();
       cell != dof_handler.end(); ++cell)
    cell->set_active_fe_index (rand() % fe_collection.size());
  
  dof_handler.distribute_dofs(fe_collection);

  std::vector<types::global_dof_index> local_dof_indices;
  for (typename hp::DoFHandler<dim>::active_cell_iterator
         cell=dof_handler.begin_active();
       cell!=dof_handler.end(); ++cell)
    {
      local_dof_indices.resize (cell->get_fe().dofs_per_cell);
      cell->get_dof_indices (local_dof_indices);

      deallog << cell << std::endl;
      for (unsigned int i=0; i<local_dof_indices.size(); ++i)
        deallog << local_dof_indices[i] << ' ';
      deallog << std::endl;
    }
}


int main ()
{
  std::ofstream logfile("random/output");
  logfile.precision(2);
  
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);  

  test<1> ();
  test<2> ();
  test<3> ();
  
  deallog << "OK" << std::endl;
}
