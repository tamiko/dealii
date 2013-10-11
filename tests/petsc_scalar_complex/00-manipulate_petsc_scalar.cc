                           

// PETSc includes
#include <petscksp.h>
#include <petscsys.h>
#include <petscmath.h>
    
#include <fstream>
#include <iostream>
#include <complex> 
#include <cassert>

// The entire point of these tests is to check that
// std::complex<double> and PetscScalar complex are compatible without
// explicitly casting. Suprisingly, it seems they are.

// Initialize a std::complex number from an PETSc complex number
void divide_petsc_complex_by_a_number ()
{
  PetscScalar alpha      = 1.0 + 2.0i;
  const PetscScalar beta = alpha/2.;

  std::cout << "   divide a petsc complex by 2.: "
	    << PetscRealPart (beta)
	    << "+" << PetscImaginaryPart (beta)
	    << "i"
	    << " should be 0.5+1i"
	    << std::endl;

  alpha /= 2.;

  std::cout << "   divide-equals a petsc complex by 2.: "
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
  const PetscScalar alpha = 1.0 + 2.0i;
  const std::complex<double> beta = alpha;

  std::cout << "   make std::complex from petsc complex: "
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
  
  std::cout << "   make petsc complex from std::complex: "
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
  
  std::cout << "   make petsc complex: "
	    << PetscRealPart (alpha)
	    << "+" << PetscImaginaryPart (alpha)
	    << "i"
	    << " should be 1+2i"
	    << std::endl;
}


// Initialize a PETSc complex number directly only.
void init_petsc_complex ()
{
  const PetscScalar alpha; 
  
  std::cout << "   make petsc complex: "
	    << PetscRealPart (alpha)
	    << "+" << PetscImaginaryPart (alpha)
	    << "i"
	    << " should be 0+0i"
	    << std::endl;

  const PetscScalar beta = 0.; 
  
  std::cout << "   make petsc complex: "
	    << PetscRealPart (beta)
	    << "+" << PetscImaginaryPart (beta)
	    << "i"
	    << " should be 0+0i"
	    << std::endl;
}               
     
int main (int argc, char **argv)
{
  try
    {
      //      PetscInitialize (&argc, &argv, (char*) 0, (char*) 0);

      {
	// some initialisations makes and so on...
	init_petsc_complex ();
	make_petsc_complex ();

	// first do conversions explicitly
	make_petsc_complex_from_std_complex ();
	make_std_complex_from_petsc_complex ();

	// then try to use operators to do the same thing
	std::cout << "   make petsc complex from std::complex: ";
	const std::complex<double> number1 (1,2);
	PetscScalar alpha1;
	alpha1 += number1; 

	std::cout << PetscRealPart (alpha1)
		  << "+" << PetscImaginaryPart (alpha1)
		  << "i"
		  << " should be 1+2i"
		  << std::endl;

	// then try to use operators to do the same thing - except it
	// seems to work without(!)
	std::cout << "   make std::complex from petsc complex: ";
	const PetscScalar alpha2 = 1.0 + 2.0i; 
	std::complex<double> number2;
	number2 += alpha2; 

	std::cout << std::real (number2)
		  << "+" << std::imag (number2)
		  << "i"
		  << " should be 1+2i"
		  << std::endl;

	std::cout << "   make std::complex from petsc complex: ";
	const PetscScalar alpha3 = 1.0 - 2.0i; 
	number2 += alpha3; 

	std::cout << std::real (number2)
		  << "+" << std::imag (number2)
		  << "i"
		  << " should be 2+0i"
		  << std::endl;

	// Let's try some other things...
	divide_petsc_complex_by_a_number ();

      }

      //      PetscFinalize ();
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

  std::cout << std::endl;

  return 0;
}


