//-----------------------------------------------------------
//
//    Copyright (C) 2017 - 2023 by the deal.II authors
//
//    This file is part of the deal.II library.
//
//    The deal.II library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE.md at
//    the top level directory of deal.II.
//
//-----------------------------------------------------------

#include <deal.II/base/parameter_handler.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/sundials/kinsol.h>

#include "../tests.h"

// Like the _05 test but in addition, this test checks that an optional
// custom_setup callback is called correctly.

namespace
{
  /**
   * Callback to test whether KINSOL::custom_setup() works.
   */
  void
  kinsol_info_callback(const char *module, const char *function, char *msg, void *ih_data)
  {
    (void)ih_data;
    deallog << "KINSOL info: " << module << ":" << function << ": " << msg << std::endl;
  }
} // namespace

int
main()
{
  initlog();

  using VectorType = Vector<double>;

  SUNDIALS::KINSOL<VectorType>::AdditionalData data;
  ParameterHandler                             prm;
  data.add_parameters(prm);

  std::ifstream ifile(SOURCE_DIR "/kinsol_linesearch.prm");
  prm.parse_input(ifile);

  // Update the Jacobian in each iteration:
  data.maximum_setup_calls = 1;


  // Size of the problem
  const unsigned int N = 2;

  SUNDIALS::KINSOL<VectorType> kinsol(data);

  kinsol.reinit_vector = [N](VectorType &v) { v.reinit(N); };

  kinsol.residual = [](const VectorType &u, VectorType &F) {
    deallog << "Evaluating the solution at u=(" << u[0] << ',' << u[1] << ')' << std::endl;

    F(0) = std::cos(u[0] + u[1]) - 1 + 2 * u[0];
    F(1) = std::sin(u[0] - u[1]) + 2 * u[1];
  };


  kinsol.iteration_function = [](const VectorType &u, VectorType &F) {
    // We want a Newton-type scheme, not a fixed point iteration. So we
    // shouldn't get into this function.
    std::abort();

    // But if anyone wanted to see how it would look like:
    F(0) = std::cos(u[0] + u[1]) - 1 + 2 * u[0] - u[0];
    F(1) = std::sin(u[0] - u[1]) + 2 * u[1] - u[1];
  };

  FullMatrix<double> J_inverse(2, 2);

  kinsol.setup_jacobian = [&J_inverse](const VectorType &u, const VectorType &F) {
    // We don't do any kind of set-up in this program, but we can at least
    // say that we're here
    deallog << "Setting up Jacobian system at u=(" << u[0] << ',' << u[1] << ')' << std::endl;

    FullMatrix<double> J(2, 2);
    J(0, 0) = -std::sin(u[0] + u[1]) + 2;
    J(0, 1) = -std::sin(u[0] + u[1]);
    J(1, 0) = std::cos(u[0] - u[1]);
    J(1, 1) = -std::cos(u[0] - u[1]) + 2;

    J_inverse.invert(J);
  };


  kinsol.solve_with_jacobian = [&J_inverse](const VectorType &rhs, VectorType &dst, double) {
    deallog << "Solving Jacobian system with rhs=(" << rhs[0] << ',' << rhs[1] << ')' << std::endl;

    J_inverse.vmult(dst, rhs);
  };

  kinsol.custom_setup = [](void *kinsol_mem) {
    // test custom_setup callback by querying some information from KINSOL
    KINSetInfoHandlerFn(kinsol_mem, kinsol_info_callback, nullptr);
    KINSetPrintLevel(kinsol_mem, 1);
  };

  VectorType v(N);
  v(0) = 0.5;
  v(1) = 1.234;

  auto niter = kinsol.solve(v);
  v.print(deallog.get_file_stream());
  deallog << "Converged in " << niter << " iterations." << std::endl;
}
