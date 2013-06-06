//----------------------------  dof_tools_8.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  dof_tools_8.cc  ---------------------------

#include "../tests.h"
#include "dof_tools_common.h"

// check
//   DoFTools::map_dof_to_boundary_indices(const DoFHandler<int>     &,
//                                         std::vector<unsigned int> &)


std::string output_file_name = "dof_tools_08/output";


template <int dim>
void
check_this (const DoFHandler<dim> &dof_handler)
{
  std::vector<types::global_dof_index> map(dof_handler.n_dofs());
  DoFTools::map_dof_to_boundary_indices (dof_handler, map);
  for (unsigned int i=0; i<map.size(); ++i)
    deallog << map[i] << " ";
  deallog << std::endl;
}
