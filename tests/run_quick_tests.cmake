## ---------------------------------------------------------------------
##
## Copyright (C) 2013 - 2016 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
##
## ---------------------------------------------------------------------

#
# This file is run when "make test" is executed by the user and is
# responsible for running the tests and printing some helpful error
# messages.
#

include(ProcessorCount)
PROCESSORCOUNT(_n_processors)

if(_n_processors EQUAL 0)
  set(_n_processors "1")
endif()

# Avoid race conditions with native Windows build tools:
if(CMAKE_HOST_SYSTEM_NAME MATCHES "Windows")
  set(_n_processors "1")
endif()

if("${CMAKE_BUILD_TYPE}" STREQUAL "")
  set(CMAKE_BUILD_TYPE "Debug")
endif()
message(STATUS "Running quick_tests in ${CMAKE_BUILD_TYPE} mode with -j${_n_processors}:")

#
# Always restrict quick tests with specified build type, but run the step
# and affinity quick tests in all available configuration:
#

string(TOLOWER "${CMAKE_BUILD_TYPE}" _build_type)
execute_process(COMMAND ${CMAKE_CTEST_COMMAND}
  -j${_n_processors} -C ${CMAKE_BUILD_TYPE} --force-new-ctest-process
  -R "quick_tests/(step.debug|step.release|affinity|.*.${_build_type})"
  OUTPUT_VARIABLE _output ERROR_VARIABLE _output RESULT_VARIABLE _return_value
  )
message(${_output})

if(NOT "${_return_value}" STREQUAL "0")
  message("
***************************************************************************
**                                                                       **
**                 Error: Some of the quick tests failed.                **
**                                                                       **
***************************************************************************

Please scroll up or check the file tests/quick_tests/quicktests.log for the
error messages. If you are unable to fix the problems, see the FAQ or write
to the mailing list linked at http://www.dealii.org\n"
    )

  string(REPLACE "\n" ";" _output "${_output}")
  foreach(_line ${_output})
    if ("${_line}" MATCHES "affinity.*FAILED")
      message("
The affinity test can fail when you are linking in a library like BLAS
which uses OpenMP. Even without calling any BLAS functions, OpenMP messes
with the thread affinity which causes TBB to run single-threaded only. You
can fix this by exporting OMP_NUM_THREADS=1. Also see GOMP_CPU_AFFINITY
and OMP_PROC_BIND.\n"
        )
    endif()

    if (${_line} MATCHES "step-petsc.*FAILED")
      message("
Additional information about PETSc issues is available
at:\nhttp://www.dealii.org/developer/external-libs/petsc.html\n"
        )
    endif()

    if (${test} MATCHES "p4est.*FAILED" AND NOT EXISTS ${test}-OK)
      message("
The p4est test can fail if you are running an OpenMPI version before 1.5.
This is a known problem and the only work around is to update to a more
recent version or use a different MPI library like MPICH.\n"
        )
    endif()
  endforeach()
endif()
