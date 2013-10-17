
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
#include <deal.II/base/logstream.h>
    
#include <fstream>
#include <iostream>
#include <cassert>

#ifdef PETSC_USE_COMPLEX

std::ofstream logfile ("00-manipulate_petsc_scalar/output");

// These tests were used while building extensions for
// PetscSclar=complex. They check that std::complex<double> and
// PetscScalar complex are compatible without explicitly casting.
//
// These tests are archaic and produce unusual verbose output for
// historical reasons...

// Divide a PETSc complex number by a real number.
void divide_petsc_complex_by_a_number ()
{
  PetscScalar alpha      = 1.0 + 2.0i;
  const PetscScalar beta = alpha/2.;

  logfile << "   divide a petsc complex by 2.: "
	  << PetscRealPart (beta)
	  << "+" << PetscImaginaryPart (beta)
	  << "i"
	  << " should be 0.5+1i"
	  << std::endl;
  
  alpha /= 2.;
  
  logfile << "   divide-equals a petsc complex by 2.: "
	  << PetscRealPart (alpha)
	  << "+" << PetscImaginaryPart (alpha)
	  << "i"
	  << " should be 0.5+1i"
	  << std::endl;
  
  assert (alpha==beta);
}

// Initialize a std::complex number from an PETSc complex number
void make_std_complex_from_petsc_complex ()
{
  const PetscScalar alpha         = 1.0 + 2.0i;
  const std::complex<double> beta = alpha;
  
  logfile << "   make std::complex from petsc complex: "
	  << std::real (beta)
	  << "+" << std::imag (beta)
	  << "i"
	  << std::endl;
  
  // even this works
  assert (alpha==beta);
}

// Initialize a PETSc complex number from an std::complex number
void make_petsc_complex_from_std_complex ()
{
  const std::complex<double> beta (1,2);
  const PetscScalar alpha = beta; 
  
  logfile << "   make petsc complex from std::complex: "
	  << PetscRealPart (alpha)
	  << "+" << PetscImaginaryPart (alpha)
	  << "i"
	  << std::endl;
  
  assert (alpha==beta);
}

// Initialize a PETSc complex number directly.
void make_petsc_complex ()
{
  const PetscScalar alpha = 1.0 + 2.0i; 
  
  logfile << "   make petsc complex: "
	  << PetscRealPart (alpha)
	  << "+" << PetscImaginaryPart (alpha)
	  << "i"
	  << " should be 1+2i"
	  << std::endl;
}


// Initialize a PETSc complex number directly only, check he is
// initialised to 0+0i.
void init_petsc_complex ()
{
  const PetscScalar alpha; 
  
  logfile << "   make petsc complex: "
	  << PetscRealPart (alpha)
	  << "+" << PetscImaginaryPart (alpha)
	  << "i"
	  << " should be 0+0i"
	  << std::endl;
  
  const PetscScalar beta = 0.; 
  
  logfile << "   make petsc complex: "
	  << PetscRealPart (beta)
	  << "+" << PetscImaginaryPart (beta)
	  << "i"
	  << " should be 0+0i"
	  << std::endl;

  assert (alpha==beta);
}               

int main (int argc, char **argv)
{
  dealii::deallog.attach (logfile);
  dealii::deallog.depth_console (1);
  
  try
    {
      PetscInitialize (&argc, &argv, (char*) 0, (char*) 0);
      {
	// some initialisations makes and so on...
	init_petsc_complex ();
	make_petsc_complex ();
	
	// first do conversions explicitly
	make_petsc_complex_from_std_complex ();
	make_std_complex_from_petsc_complex ();
	
	// then try to use operators to do the same thing
	logfile << "   make petsc complex from std::complex: ";
	const std::complex<double> number1 (1,2);
	PetscScalar alpha1;
	alpha1 += number1; 

	logfile << PetscRealPart (alpha1)
		<< "+" << PetscImaginaryPart (alpha1)
		<< "i"
		<< " should be 1+2i"
		<< std::endl;
	
	// then try to use operators to do the same thing - except it
	// seems to work without(!)
	logfile << "   make std::complex from petsc complex: ";
	const PetscScalar alpha2 = 1.0 + 2.0i; 
	std::complex<double> number2;
	number2 += alpha2; 
	
	logfile << std::real (number2)
		<< "+" << std::imag (number2)
		<< "i"
		<< " should be 1+2i"
		<< std::endl;

	logfile << "   make std::complex from petsc complex: ";
	const PetscScalar alpha3 = 1.0 - 2.0i; 
	number2 += alpha3; 
	
	logfile << std::real (number2)
		<< "+" << std::imag (number2)
		  << "i"
		<< " should be 2+0i"
		<< std::endl;
	
	// Let's try some other things...
	divide_petsc_complex_by_a_number ();

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

#endif // PETSC_USE_COMPLEX
