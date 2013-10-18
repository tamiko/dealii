
// PETSc includes: We would like to do this cleanly
//
#include <petscksp.h>
#include <petscsys.h>
#include <petscmath.h>
//
// though a horrible warning message 
//
// warning: imaginary constants are a GCC extension [enabled by default]
//
// appears. This seems to be related to the direct useage of
// PetscScalar in accordance with gcc.gnu.org bug 7263.

// deal.II includes
#include "../tests.h"
#include <deal.II/base/logstream.h>
    
#include <fstream>
#include <iostream>
#include <cassert>

// These tests were used while building extensions for
// PetscSclar=complex. They check that std::complex<double> and
// PetscScalar complex are compatible without explicitly casting.
//
// These tests look archaic - they are the original tests used by the
// authors to figure out how PetscScalar complex works...

// Multiply a PETSc complex number by a complex number.
void multiply_petsc_complex_by_a_number ()
{
  deallog << "Check PetscScalar multiply" << std::endl;

  // make alpha
  const PetscScalar alpha = 1.0 + 2.0*PETSC_i;

  // make beta 
  const PetscScalar beta  = 2.0 + 1.0*PETSC_i;
  
  // multiply them together
  const PetscScalar gamma = alpha * beta;

  // Check real access
  Assert (PetscRealPart (gamma)==0., ExcInternalError());
  Assert (PetscImaginaryPart (gamma) ==5., ExcInternalError());

  // This is legal too!
  Assert (gamma==std::complex<double> (0.,5), ExcInternalError());

  deallog << "OK" << std::endl;
}


// Divide a PETSc complex number by a real number.
void divide_petsc_complex_by_a_number ()
{
  deallog << "Check PetscScalar division" << std::endl;

  // make alpha
  PetscScalar alpha      = 1.0 + 2.0*PETSC_i;

  // make beta <- alpha/2
  const PetscScalar beta = alpha/2.;

  // Check real access
  Assert (PetscRealPart (alpha)==1.,  ExcInternalError());
  Assert (PetscRealPart (beta) ==0.5, ExcInternalError());
  
  // operate on alpha (now alpha==beta)
  alpha /= 2.;
  
  // Check complex access
  Assert (PetscImaginaryPart (alpha)==1., ExcInternalError());
  Assert (PetscImaginaryPart (beta) ==1., ExcInternalError());

  // This is legal too!
  Assert (alpha==beta, ExcInternalError());

  deallog << "OK" << std::endl;
}

// Initialize a std::complex number from an PETSc complex number
void make_std_complex_from_petsc_complex ()
{
  deallog << "Check std initialised from PetscScalar" << std::endl;
  const PetscScalar alpha         = 1.0 + 2.0*PETSC_i;
  const std::complex<double> beta = alpha;
  
  // These should be the same of course
  Assert (alpha==beta, ExcInternalError());

  deallog << "OK" << std::endl;
}

// Initialize a PETSc complex number from an std::complex number
void make_petsc_complex_from_std_complex ()
{
  deallog << "Check PetscScalar initialised from std" << std::endl;
  const std::complex<double> beta (1,2);
  const PetscScalar alpha = beta; 

  // These should be the same of course
  Assert (alpha==beta, ExcInternalError());

  deallog << "OK" << std::endl;
}

// Initialize a PETSc complex number directly.
void make_petsc_complex ()
{
  deallog << "Check a nonzero PetscScalar" << std::endl;

  // init
  const PetscScalar alpha = 1.0 + 2.0*PETSC_i; 

  // Test if PetscScalar is an std::complex.
  Assert (alpha==std::complex<double> (1.,2.), ExcInternalError());

  deallog << "OK" << std::endl;
}


// Initialize a PETSc complex number directly only, check he is
// initialised to 0+0i.
void init_petsc_complex ()
{
  deallog << "Check PetscScalar initialisation" << std::endl;

  // Initialise (no argument) to zero.
  const PetscScalar alpha; 
  Assert (alpha==std::complex<double> (0.,0.), ExcInternalError());
  
  // Initialise (real argument) to zero.
  const PetscScalar beta = 0.; 
  Assert (beta==alpha, ExcInternalError());

  // Initialise (real+complex argument) to zero.
  const PetscScalar gamma = 0. + 0.*PETSC_i; 
  Assert (gamma==beta, ExcInternalError());

  // If alpha==beta==gamma, then:
  Assert (alpha==gamma, ExcInternalError());

  deallog << "OK" << std::endl;
}               

int main (int argc, char **argv)
{
  std::ofstream logfile ("00/output");
  dealii::deallog.attach (logfile);
  dealii::deallog.depth_console (0);
  deallog.threshold_double(1.e-10);
  
  try
    {
      PetscInitialize (&argc, &argv, (char*) 0, (char*) 0);
      {
	// check initialisation (zero and nonzero)
	init_petsc_complex ();
	make_petsc_complex ();
	
	// check initialisation from std::complex (and vice versa)
	make_petsc_complex_from_std_complex ();
	make_std_complex_from_petsc_complex ();
		
	// Check some operators
	divide_petsc_complex_by_a_number ();
	multiply_petsc_complex_by_a_number ();
      }
      PetscFinalize ();
    }
  
  
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  
  logfile << std::endl;
  
  return 0;
}


