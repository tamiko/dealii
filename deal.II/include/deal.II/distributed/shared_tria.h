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
  namespace shared
  {

   /**
    * This is an extension of dealii::Triangulation class to automatically
    * partition trianguation using MPI.
    */
    template <int dim, int spacedim = dim>
    class Triangulation : public dealii::Triangulation<dim,spacedim>
    {
    public:
      typedef typename dealii::Triangulation<dim,spacedim>::active_cell_iterator active_cell_iterator;
      typedef typename dealii::Triangulation<dim,spacedim>::cell_iterator        cell_iterator;

      /**
       * Constructor.
       */
      Triangulation (MPI_Comm mpi_communicator);

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
      
      
    private:
      MPI_Comm mpi_communicator;
      types::subdomain_id my_subdomain;
      types::subdomain_id num_subdomains;
    };
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
