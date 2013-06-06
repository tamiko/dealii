//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2007, 2008, 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------


// there was a bug in getting the divergence of shape functions for the
// SymmetricTensor extractors. test that it is fixed by comparing with
// get_function_divergences

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <fstream>



template<int dim>
void test (const Triangulation<dim>& tr,
	   const FiniteElement<dim>& fe)
{
  deallog << "FE=" << fe.get_name()
	  << std::endl;

  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);

  Vector<double> fe_function(dof.n_dofs());
  for (unsigned int i=0; i<dof.n_dofs(); ++i)
    fe_function(i) = i+1;
  
  const QGauss<dim> quadrature(2);
  FEValues<dim> fe_values (fe, quadrature,
			   update_values | update_gradients);
  fe_values.reinit (dof.begin_active());

				   // let the FEValues object compute the
				   // divergences at quadrature points
  std::vector<Tensor<1,dim> > divergences (quadrature.size());
  FEValuesExtractors::SymmetricTensor<2> extractor(0);
  fe_values[extractor]
    .get_function_divergences (fe_function, divergences);

				   // now do the same "by hand"
  std::vector<types::global_dof_index> local_dof_indices (fe.dofs_per_cell);
  dof.begin_active()->get_dof_indices (local_dof_indices);
  
  for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
    {
      deallog << "i=" << i << std::endl;
      
      for (unsigned int q=0; q<quadrature.size(); ++q)
	deallog << "  q_point=" << q << std::endl
		<< "    div= " << fe_values[extractor].divergence (i,q) << std::endl;
    }
}



template<int dim>
void test_hyper_sphere()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_ball(tr);

  static const HyperBallBoundary<dim> boundary;
  tr.set_boundary (0, boundary);

  FESystem<dim> fe (FE_Q<dim>(1),
		    SymmetricTensor<2,dim>::n_independent_components);
  test(tr, fe);
}


int main()
{
  std::ofstream logfile ("fe_values_view_24/output");
  deallog << std::setprecision (3);

  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-7);

  test_hyper_sphere<2>();
  test_hyper_sphere<3>();
}
