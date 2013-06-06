//----------------------------  renumber_block_wise_02.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2003, 2004, 2009, 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  renumber_block_wise_02.cc  ---------------------------


// Test DoFRenumbering::block_wise. For the element used here, it
// needs to produce the exact same numbering as that for
// DoFRenumber::component_wise



#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/multigrid/mg_dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/hp/fe_collection.h>

#include <fstream>



template <int dim>
std::vector<types::global_dof_index>
get_dofs (const hp::DoFHandler<dim> &dof)
{
  std::vector<types::global_dof_index> local;
  std::vector<types::global_dof_index> global;
  for (typename hp::DoFHandler<dim>::active_cell_iterator cell=dof.begin_active();
       cell != dof.end(); ++cell)
    {
      local.resize (cell->get_fe().dofs_per_cell);
      cell->get_dof_indices (local);

      global.insert (global.end(), local.begin(), local.end());
    }

  return global;
}



template <int dim>
void
check_renumbering(hp::DoFHandler<dim>& dof)
{
				   // Prepare a reordering of
				   // components so that each
				   // component maps to its natural
				   // block
  std::vector<unsigned int> order(dof.get_fe().n_components());
  order[0] = 0;
  order[1] = 1;
  order[2] = 1;

				   // do component-wise and save the
				   // results
  DoFRenumbering::component_wise (dof, order);
  const std::vector<types::global_dof_index> vc = get_dofs (dof);

				   // now do the same with blocks
  DoFRenumbering::block_wise (dof);
  const std::vector<types::global_dof_index> vb = get_dofs (dof);

  Assert (vc == vb, ExcInternalError());

  deallog << "OK" << std::endl;
}


template <int dim>
void
check ()
{
  Triangulation<dim> tr;
  if (dim==2)
    GridGenerator::hyper_ball(tr, Point<dim>(), 1);
  else
    GridGenerator::hyper_cube(tr, -1,1);
  tr.refine_global (1);
  tr.begin_active()->set_refine_flag ();
  tr.execute_coarsening_and_refinement ();
  if (dim==1)
    tr.refine_global(2);

  hp::DoFHandler<dim> dof(tr);
  {
    bool coin = false;
    for (typename hp::DoFHandler<dim>::active_cell_iterator cell=dof.begin_active();
	 cell != dof.end(); ++cell)
      {
	cell->set_active_fe_index (coin ? 0 : 1);
	coin = !coin;
      }
  }

				   // note that the following elements
				   // have 3 components but 2 blocks
  FESystem<dim> e1 (FE_DGQ<dim>(1), 1, FESystem<dim>(FE_DGQ<dim>(2), 2), 1);
  FESystem<dim> e2 (FE_DGQ<dim>(2), 1, FESystem<dim>(FE_DGQ<dim>(1), 2), 1);

  hp::FECollection<dim> fe_collection;
  fe_collection.push_back (e1);
  fe_collection.push_back (e2);

  dof.distribute_dofs(fe_collection);
  check_renumbering(dof);
  dof.clear();
}


int main ()
{
  std::ofstream logfile ("renumber_block_wise_02/output");
  deallog << std::setprecision (2);
  deallog << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console (0);

  deallog.push ("2d");
  check<2> ();
  deallog.pop ();
  deallog.push ("3d");
  check<3> ();
  deallog.pop ();
}
