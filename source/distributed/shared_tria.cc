// ---------------------------------------------------------------------
// $Id: tria.cc 32807 2014-04-22 15:01:57Z heister $
//
// Copyright (C) 2008 - 2013 by the deal.II authors
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


#include <deal.II/base/utilities.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/distributed/tria.h>


#include <algorithm>
#include <numeric>
#include <iostream>
#include <fstream>


DEAL_II_NAMESPACE_OPEN


namespace parallel
{
  namespace shared
  {

	template <int dim, int spacedim>
	Triangulation<dim,spacedim>::NumberCache::NumberCache()
	  :
	  n_global_active_cells(0),
	  n_global_levels(0)
	 {}

    template <int dim, int spacedim>
    Triangulation<dim,spacedim>::Triangulation (MPI_Comm mpi_communicator,
    		                                    const typename dealii::Triangulation<dim,spacedim>::MeshSmoothing smooth_grid):
        dealii::parallel::Triangulation<dim,spacedim>(smooth_grid,false),
    	mpi_communicator (Utilities::MPI::
                          duplicate_communicator(mpi_communicator)),
        my_subdomain (Utilities::MPI::this_mpi_process (this->mpi_communicator)),
        num_subdomains(Utilities::MPI::n_mpi_processes(mpi_communicator))
    {
    	number_cache.n_locally_owned_active_cells.resize (num_subdomains);
    }


    template <int dim, int spacedim>
    Triangulation<dim,spacedim>::~Triangulation ()
    {

    }

    template <int dim, int spacedim>
    unsigned int
    Triangulation<dim,spacedim>::n_locally_owned_active_cells () const
    {
    	return number_cache.n_locally_owned_active_cells[my_subdomain];
    }

    template <int dim, int spacedim>
    unsigned int
    Triangulation<dim,spacedim>::n_global_levels () const
    {
    	return number_cache.n_global_levels;
    }

    template <int dim, int spacedim>
    types::global_dof_index
    Triangulation<dim,spacedim>::n_global_active_cells () const
    {
    	return number_cache.n_global_active_cells;
    }

    template <int dim, int spacedim>
    const std::vector<unsigned int> &
    Triangulation<dim,spacedim>::n_locally_owned_active_cells_per_processor () const
    {
    	return number_cache.n_locally_owned_active_cells;
    }

    template <int dim, int spacedim>
    void
    Triangulation<dim,spacedim>::update_number_cache ()
    {
    	Assert (number_cache.n_locally_owned_active_cells.size()
    			==
				Utilities::MPI::n_mpi_processes (mpi_communicator),
				ExcInternalError());

    	std::fill (number_cache.n_locally_owned_active_cells.begin(),
    			number_cache.n_locally_owned_active_cells.end(),
    			0);

    	if (this->n_levels() > 0)
    		for (typename Triangulation<dim,spacedim>::active_cell_iterator
    				cell = this->begin_active();
    				cell != this->end(); ++cell)
    			if (cell->subdomain_id() == my_subdomain)
    				++number_cache.n_locally_owned_active_cells[my_subdomain];

    	unsigned int send_value
    	= number_cache.n_locally_owned_active_cells[my_subdomain];
    	MPI_Allgather (&send_value,
    			1,
    			MPI_UNSIGNED,
    			&number_cache.n_locally_owned_active_cells[0],
    			1,
    			MPI_UNSIGNED,
    			mpi_communicator);

    	number_cache.n_global_active_cells
    	= std::accumulate (number_cache.n_locally_owned_active_cells.begin(),
    			number_cache.n_locally_owned_active_cells.end(),
    			/* ensure sum is computed with correct data type:*/
    			static_cast<types::global_dof_index>(0));
    	number_cache.n_global_levels = Utilities::MPI::max(this->n_levels(), mpi_communicator);
    }



    template <int dim, int spacedim>
    types::subdomain_id
    Triangulation<dim,spacedim>::locally_owned_subdomain () const
    {
      return my_subdomain;
    }
    
    template <int dim, int spacedim>
    void 
    Triangulation<dim,spacedim>::execute_coarsening_and_refinement () {
    	  dealii::Triangulation<dim,spacedim>::execute_coarsening_and_refinement ();
    	  dealii::GridTools::partition_triangulation (num_subdomains, *this);
    	  update_number_cache ();
    }
    
    template <int dim, int spacedim>
    void 	
    Triangulation<dim,spacedim>::create_triangulation (const std::vector< Point< spacedim > > &vertices, 
    												   const std::vector< CellData< dim > > &cells, 
    												   const SubCellData &subcelldata) {												   
      try
        {
          dealii::Triangulation<dim,spacedim>::
          create_triangulation (vertices, cells, subcelldata);
        }
      catch (const typename dealii::Triangulation<dim,spacedim>::DistortedCellList &)
        {
          // the underlying triangulation should not be checking for distorted
          // cells
          AssertThrow (false, ExcInternalError());
        }
      dealii::GridTools::partition_triangulation (num_subdomains, *this);
      update_number_cache ();
    }


    template <int dim, int spacedim>
    MPI_Comm
    Triangulation<dim,spacedim>::get_communicator () const
    {
      return mpi_communicator;
    }

  }
}


/*-------------- Explicit Instantiations -------------------------------*/
#include "shared_tria.inst"


DEAL_II_NAMESPACE_CLOSE
