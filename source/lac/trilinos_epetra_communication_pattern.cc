// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#include <deal.II/lac/trilinos_epetra_communication_pattern.h>

#ifdef DEAL_II_WITH_TRILINOS

#  include <deal.II/base/index_set.h>

#  include <Epetra_Map.h>

#  include <memory>

DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  namespace EpetraWrappers
  {
    CommunicationPattern::CommunicationPattern(const IndexSet &locally_owned_indices,
                                               const IndexSet &ghost_indices,
                                               const MPI_Comm  communicator)
    {
      // virtual functions called in constructors and destructors never use the
      // override in a derived class
      // for clarity be explicit on which function is called
      CommunicationPattern::reinit(locally_owned_indices, ghost_indices, communicator);
    }



    void
    CommunicationPattern::reinit(const IndexSet &locally_owned_indices,
                                 const IndexSet &ghost_indices,
                                 const MPI_Comm  communicator)
    {
      comm = std::make_shared<const MPI_Comm>(communicator);

      Epetra_Map vector_space_vector_map = locally_owned_indices.make_trilinos_map(*comm, false);
      Epetra_Map read_write_vector_map   = ghost_indices.make_trilinos_map(*comm, true);

      // Target map is read_write_vector_map
      // Source map is vector_space_vector_map. This map must have uniquely
      // owned GID.
      importer = std::make_unique<Epetra_Import>(read_write_vector_map, vector_space_vector_map);
    }



    MPI_Comm
    CommunicationPattern::get_mpi_communicator() const
    {
      return *comm;
    }



    const Epetra_Import &
    CommunicationPattern::get_epetra_import() const
    {
      return *importer;
    }
  } // namespace EpetraWrappers
} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#endif
