// ---------------------------------------------------------------------
// $Id: tria.h 32739 2014-04-08 16:39:47Z denis.davydov $
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

#ifndef __deal2__distributed__shared_tria_h
#define __deal2__distributed__shared_tria_h


#include <deal.II/base/config.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/grid/tria.h>

#include <deal.II/base/std_cxx1x/function.h>
#include <deal.II/base/std_cxx1x/tuple.h>

#include <set>
#include <vector>
#include <list>
#include <utility>

#ifdef DEAL_II_WITH_MPI
#  include <mpi.h>
#endif


DEAL_II_NAMESPACE_OPEN

template <int, int> class Triangulation;


namespace parallel
{
    template <int dim, int spacedim = dim>
    class Triangulation : public dealii::Triangulation<dim,spacedim>
    {
    public:

      /**
       * Constructor.
       */
      Triangulation (const typename dealii::Triangulation<dim,spacedim>::MeshSmoothing smooth_grid = dealii::Triangulation<dim,spacedim>::none,
    		         const bool check_for_distorted_cells = false) :
    		  dealii::Triangulation<dim,spacedim>(smooth_grid,check_for_distorted_cells) {};

      /**
       * Destructor.
       */
      virtual ~Triangulation () {};

      /**
       * Return locally owned subdomain id,
       * which is equivalent to the rank of the current mpi process in
       * communicator provided to class constructor.
       */
      virtual types::subdomain_id locally_owned_subdomain () const = 0;

      /**
       * Return MPI communicator used by this triangulation.
       */
      virtual MPI_Comm get_communicator () const = 0;

      /**
       * Coarsen and refine the mesh according to refinement and
       * coarsening flags set.
       *
       * This step is equivalent to the dealii::Triangulation class
       * with an addition of calling dealii::GridTools::partition_triangulation() at the end.
       */
      virtual void execute_coarsening_and_refinement () = 0;

      /**
        * Create a triangulation.
        *
        * This function also partitions triangulation based on the
        * MPI communicator provided to constructor.
        */
      virtual void 	create_triangulation (const std::vector< Point< spacedim > > &vertices,
      									  const std::vector< CellData< dim > > &cells,
      									  const SubCellData &subcelldata) = 0;

      /**
       * Return the number of active cells owned by each of the MPI
       * processes that contribute to this triangulation. The element
       * of this vector indexed by locally_owned_subdomain() equals
       * the result of n_locally_owned_active_cells().
       */
      virtual const std::vector<unsigned int> &
      n_locally_owned_active_cells_per_processor () const = 0;

      /**
       * Return the number of active cells in the triangulation that
       * are locally owned, i.e. that have a subdomain_id equal to
       * locally_owned_subdomain().
       *
       */
      virtual unsigned int n_locally_owned_active_cells () const = 0;

      /**
       * Return the sum over all processors of the number of active
       * cells owned by each processor. This equals the overall number
       * of active cells in the shared triangulation.
       */
      virtual types::global_dof_index n_global_active_cells () const = 0;

      /**
       * Returns the global maximum level. This may be bigger than n_levels.
       */
      virtual unsigned int n_global_levels () const = 0;
    };


  namespace shared
  {

   /**
    * This is an extension of dealii::Triangulation class to automatically
    * partition trianguation using MPI.
    */
    template <int dim, int spacedim = dim>
    class Triangulation : public dealii::parallel::Triangulation<dim,spacedim>
    {
    public:
      typedef typename dealii::Triangulation<dim,spacedim>::active_cell_iterator active_cell_iterator;
      typedef typename dealii::Triangulation<dim,spacedim>::cell_iterator        cell_iterator;

      /**
       * Constructor.
       */
      Triangulation (MPI_Comm mpi_communicator,
    		         const typename dealii::Triangulation<dim,spacedim>::MeshSmoothing =
		     (dealii::Triangulation<dim,spacedim>::none) );

      /**
       * Destructor.
       */
      virtual ~Triangulation ();

      /**
       * Return locally owned subdomain id,
       * which is equivalent to the rank of the current mpi process in
       * communicator provided to class constructor.
       */
      types::subdomain_id locally_owned_subdomain () const;

      /**
       * Return MPI communicator used by this triangulation.
       */
      MPI_Comm get_communicator () const;
      
      /**
       * Coarsen and refine the mesh according to refinement and
       * coarsening flags set.
       *
       * This step is equivalent to the dealii::Triangulation class
       * with an addition of calling dealii::GridTools::partition_triangulation() at the end.
       */
      virtual void execute_coarsening_and_refinement ();
      
      /**
        * Create a triangulation.
        *
        * This function also partitions triangulation based on the
        * MPI communicator provided to constructor.
        */
      virtual void 	create_triangulation (const std::vector< Point< spacedim > > &vertices, 
      									  const std::vector< CellData< dim > > &cells, 
      									  const SubCellData &subcelldata);
      
      /**
       * Return the number of active cells owned by each of the MPI
       * processes that contribute to this triangulation. The element
       * of this vector indexed by locally_owned_subdomain() equals
       * the result of n_locally_owned_active_cells().
       */
      const std::vector<unsigned int> &
      n_locally_owned_active_cells_per_processor () const;

      /**
       * Return the number of active cells in the triangulation that
       * are locally owned, i.e. that have a subdomain_id equal to
       * locally_owned_subdomain().
       *
       */
      unsigned int n_locally_owned_active_cells () const;

      /**
       * Return the sum over all processors of the number of active
       * cells owned by each processor. This equals the overall number
       * of active cells in the shared triangulation.
       */
      types::global_dof_index n_global_active_cells () const;

      /**
       * Returns the global maximum level. This may be bigger than n_levels.
       */
      virtual unsigned int n_global_levels () const;

      
    private:
      MPI_Comm mpi_communicator;
      types::subdomain_id my_subdomain;
      types::subdomain_id num_subdomains;

      struct NumberCache
      {
    	  std::vector<unsigned int> n_locally_owned_active_cells;
    	  types::global_dof_index   n_global_active_cells;
    	  unsigned int              n_global_levels;

    	  NumberCache();
      };

      NumberCache number_cache;

      /**
       * Update the number_cache variable after mesh creation or
       * refinement.
       */
      void update_number_cache ();


    };
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
