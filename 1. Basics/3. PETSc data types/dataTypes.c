static char help[] = "Hello, world!\n\n";

#include <petscsys.h>

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
    /*
       PetscInt is a PETSc type that represents an integer. Used primarily to
       represent size of arrays and indexing into arrays. Its size can be
       configured with the option --with-64-bit-indices to be either 32bit or
       64bit.
    */
    PetscInt integer = 1;

    /*
       PetscScalar is a PETSc type that represents either a double precision
       number, a double precision complex number, a single precision real
       number, a long double, or an int.
    */
    PetscScalar scalar = 1.0;

    /*
       PetscReal is a PETSc type that represents a real number version of
       PetscScalar.
    */
    PetscReal real = 1.0;

    /*
       PetscComplex is a PETSc type that represents a complex number with
       precision matching that of PetscReal. Complex numbers are automatically
       available if PETSc was able to find a working complex implementation.
    */
    PetscComplex complexNumber = 1.0 + 2.0 * PETSC_i;

    /*
       PetscBool is a logical variable. It is actually an int in C and a
       logical in Fortran.
    */
    PetscBool bool ;

    PetscErrorCode ierr;

    PetscInitialize(&argc, &argv, (char*) 0, help);

    PetscFinalize();

    return 0;
}
