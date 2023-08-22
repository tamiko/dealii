// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2021 by the deal.II authors
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


// Check that one can output physical units along with output
// fields. This test is for vector-valued quantities and checks the
// pvtu writer.


#include <deal.II/base/data_out_base.h>

#include <string>
#include <vector>

#include "../tests.h"

#include "patches.h"

// Output data on repetitions of the unit hypercube

template <int dim, int spacedim>
void
check(DataOutBase::VtkFlags flags, std::ostream &out)
{
  std::vector<std::string> names(5);
  names[0] = "x1";
  names[1] = "x2";
  names[2] = "x3";
  names[3] = "x4";
  names[4] = "i";
  for (unsigned int i = 1; i < 1 + dim; ++i)
    names[i] = "v";

  using Descriptor =
    std::tuple<unsigned int, unsigned int, std::string, DataComponentInterpretation::DataComponentInterpretation>;
  std::vector<Descriptor> vectors(
    1, Descriptor{1, 1 + dim - 1, "v", DataComponentInterpretation::component_is_part_of_vector});
  DataOutBase::write_pvtu_record(out, {"abc.vtu"}, names, vectors, flags);
}


template <int dim, int spacedim>
void
check_all(std::ostream &log)
{
  DataOutBase::VtkFlags flags;

  // Set the units of some but not all output quantities
  flags.physical_units["x1"] = "kg/s";
  flags.physical_units["v"]  = "m/s";

  log << "==============================" << std::endl
      << dim << spacedim << ".vtu" << std::endl
      << "==============================" << std::endl;
  check<dim, spacedim>(flags, log);
}

int
main()
{
  std::stringstream ss;
  check_all<1, 1>(ss);
  check_all<1, 2>(ss);
  check_all<2, 2>(ss);
  check_all<2, 3>(ss);
  check_all<3, 3>(ss);

  std::ofstream logfile("output");
  filter_out_xml_key(ss, "DataArray", logfile);
}
