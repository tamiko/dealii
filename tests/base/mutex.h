// ---------------------------------------------------------------------
//
// Copyright (C) 2023 by the deal.II authors
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

// Check that our std::mutex wrapper works as intended


#include <deal.II/base/mutex.h>

#include "../tests.h"

int
main()
{
  initlog();

  Threads::Mutex mutex;
  mutex.lock();

  // copy constructor
  Threads::Mutex mutex2(mutex);
  mutex2.lock();
  mutex2.unlock();

  // copy assignment
  mutex2 = mutex;
  mutex2.lock();
  mutex2.unlock();

  mutex.unlock();

  // move constructor
  Threads::Mutex mutex3(Threads::Mutex{});
  mutex3.lock();
  mutex3.unlock();

  // move assignment
  mutex3 = Threads::Mutex{};
  mutex3.lock();
  mutex3.unlock();

  deallog << "OK!" << std::endl;
}

