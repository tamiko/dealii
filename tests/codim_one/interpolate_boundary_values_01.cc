
//----------------------------  template.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005, 2008, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  template.cc  ---------------------------


// test VectorTools::interpolate_boundary_values for codim=1

#include "../tests.h"
#include <fstream>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/vector_tools.h>

#include <string>

std::ofstream logfile("interpolate_boundary_values_01/output");

template <int dim, int spacedim>
void test(std::string filename) {
    Triangulation<dim, spacedim> tria;
    GridIn<dim, spacedim> gi;
    gi.attach_triangulation (tria);
    std::ifstream in (filename.c_str());
    gi.read_ucd (in);

    deallog << tria.n_active_cells() << " active cells" << std::endl;

    FE_Q<dim,spacedim> fe(2);
    DoFHandler<dim,spacedim> dof_handler (tria);
    dof_handler.distribute_dofs (fe);

    deallog << dof_handler.n_dofs() << " degrees of freedom" << std::endl;

    std::map<types::global_dof_index, double> bv;
    VectorTools::interpolate_boundary_values (dof_handler,
					      0,
					      Functions::SquareFunction<spacedim>(),
					      bv);
    deallog << bv.size() << " boundary degrees of freedom" << std::endl;

    for (std::map<types::global_dof_index, double>::const_iterator i = bv.begin();
	 i != bv.end(); ++i)
      deallog << i->first << ' ' << i->second << std::endl;

    for (typename DoFHandler<dim,spacedim>::active_cell_iterator
	   cell = dof_handler.begin_active(); cell != dof_handler.end(); ++cell)
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
	if (cell->at_boundary(f))
	  for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_face; ++v)
	    for (unsigned int i=0; i<fe.dofs_per_vertex; ++i)
	      {
		Assert (bv.find(cell->face(f)->vertex_dof_index(v,i))
			!= bv.end(),
			ExcInternalError());
		Assert (bv[cell->face(f)->vertex_dof_index(v,i)]
			==
			Functions::SquareFunction<spacedim>()
			.value(cell->face(f)->vertex(v),i),
			ExcInternalError());
	      }
}



int main ()
{
  deallog.attach(logfile);
  deallog.depth_console(0);

  test<2,3>("grids/square.inp");
  test<2,3>("grids/sphere_1.inp");

  return 0;
}

